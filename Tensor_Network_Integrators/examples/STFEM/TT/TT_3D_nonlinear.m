close all; clear; clc;

addpath(genpath('../../../matlab/STFEM/src/'))
addpath(genpath('../../../matlab/utils/chebfun/'))
addpath(genpath('../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../matlab/utils/ttfunc/'))

%% problem definitions
% du/dt + DiffU = u - u.^3

uexactfn = @(t,x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z).*sin(pi*t) ...
  + sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z).*sin(2*pi*t);

afn = @(x,y,z) 1.0;

% compute terms
dudt = compute_dudt_fn(uexactfn);
Dfn = compute_3D_diffusion_fn(uexactfn,afn);
gfn = @(t,x,y,z) dudt(t,x,y,z) + Dfn(t,x,y,z)...
- uexactfn(t,x,y,z) + uexactfn(t,x,y,z).^3;
%%
d = 3;
NEs = 1 + 2.^[3:4]; %full run [3:9]

R = cell(numel(NEs),1);
fname='../plot_data/TT_3D_nonlinear.mat';

for jjns = 1:numel(NEs)
  %% parameters
  ts = datetime;
  Lx = 1; % Length of the domain in x-direction
  Ly = 1; % Length of the domain in y-direction
  Lz = 1;
  Lt = 1;

  Ex = NEs(jjns); % Number of elements in x-direction
  hx = Lx / Ex; % Element size in x-direction
  Nx = Ex + 1; % Number of nodes in x-direction

  Et = Ex-1;
  ht = Lt/Et;
  Nt = Et + 1;

  tt_tol = 0.01*hx^2;
  % tt_tol = 1e-8;
  %% Create mesh
  X = 0:hx:1;
  T = 0:ht:1;

  %% create grid in tt
  Itt = [ones(1,Nt),repmat({ones(1, Nx)}, 1, d)];
  Ctt = cell(1, d+1);
  temp = Itt;
  temp{1} = T;
  Ctt{1} = cell2core(tt_tensor,temp);

  for ic = 2:d+1
    temp = Itt;
    temp{ic} = X;
    Ctt{ic} = cell2core(tt_tensor,temp);
  end

  %% create spatial grid in tt
  Itt = repmat({ones(1, Nx)}, 1, d);
  Cxtt = cell(1, d);
  for ic = 1:d
    temp = Itt;
    temp{ic} = X;
    Cxtt{ic} = cell2core(tt_tensor,temp);
  end
  %% %%%%%%%%%%% Construction of the left hand side matrix
  % Initialize NNL with zeros
  NNL = sparse(Ex, 2*Ex);  % for space dimension
  % Populate NNL matrix
  for ii=2:Ex
    NNL(ii,ii+(ii-2))=1; NNL(ii,ii+(ii-2)+1)=1;
  end
  NNL(1,1)=1; NNL(Nx, 2*Ex)=1;

  NNR=NNL';

  NNLt = sparse(Et, 2*Et);  % for time dimension
  % Populate NNLt matrix
  for ii=2:Et
    NNLt(ii,ii+(ii-2))=1; NNLt(ii,ii+(ii-2)+1)=1;
  end
  NNLt(1,1)=1; NNLt(Nt, 2*Et)=1;

  NNRt=NNLt';


  % Assemble stiffness matrix and load vector
  [A1,A2,M1,M2,B1,B2] = nonlinear_Mat(X(1:2));


  AA = {kron(eye(Ex),A1),kron(eye(Ex),A2)};
  MM = {kron(eye(Ex),M1),kron(eye(Ex),M2)};
  BB = {kron(eye(Ex),B1),kron(eye(Ex),B2)};
  % Element stiffness matrix for a Q1 element in 2D
  Bs=(hx)*(1/6)*[2 1;1 2];
  ML = NNL*kron(eye(Ex),Bs)*NNR;


  % matrix for time dimension
  [A1t,A2t,M1t,M2t,B1t,B2t] = nonlinear_Mat(T(1:2));
  Bt=(ht)*(1/6)*[2 1;1 2];
  MT = Mat_Time(T(1:2));
  ZT = NNLt*kron(eye(Et),MT)*NNRt;
  MLt = NNLt*kron(eye(Et),Bt)*NNRt;
  MX1t = NNLt*kron(eye(Et),M1t)*NNRt;
  MX2t = NNLt*kron(eye(Et),M2t)*NNRt;

  %% Build global operator for time
  Ix = [2:Ex]; % interior system
  It = [2:Et+1];
  AT = matrices_to_tt_matrix_fn({ZT(It,It),ML(Ix,Ix),ML(Ix,Ix),ML(Ix,Ix)});
  %%  Build Laplace Operator
  % nterm = 2; %number of rank-1 term in the decomposition of funciton a
  att= amen_cross_zero(Cxtt, @(x) cross_fun_nD(x,afn),tt_tol,'verb',0);
  %
  G = core2cell(att);
  nt1 = att.r(2);
  nt2 = att.r(3);

  Agtt = [];
  for j1 = 1:nt1
    for j2 = 1:nt2
      % Build 1D matrix operator
      CC = {kron(G{1}(1,:,j1)', [1;1]), kron(G{2}(j1,:,j2)', [1;1]), kron(G{3}(j2,:,1)', [1;1])};

      % Calculate diagonal matrices only once
      diag_CC = cell(d,2);
      for idim = 1:d
        diag_CC{idim,1} = diag(CC{idim}(1:end-2));
        diag_CC{idim,2} = diag(CC{idim}(3:end));
      end

      % Compute AG and MG matrices with diag multiplications
      for idim = 1:d
        for ipt = 1:2
          AG{idim,ipt} = NNL * (AA{ipt} * diag_CC{idim,ipt}) * NNR;
          MG{idim,ipt} = NNL * (MM{ipt} * diag_CC{idim,ipt}) * NNR;
        end
      end

      % Calculate the current TT matrix
      Agttcur =[];
      for ipt = 1:2
        for ipt2 = 1:2
          for ipt3 = 1:2
            curterm = round(...
              matrices_to_tt_matrix_fn({AG{1,ipt}(Ix, Ix), ...
              MG{2,ipt2}(Ix, Ix), MG{3,ipt3}(Ix, Ix)}) ...
              + matrices_to_tt_matrix_fn({MG{1,ipt}(Ix, Ix), ...
              AG{2,ipt2}(Ix, Ix),MG{3,ipt3}(Ix, Ix)}) ...
              + matrices_to_tt_matrix_fn({MG{1,ipt}(Ix, Ix), ...
              MG{2,ipt2}(Ix, Ix),AG{3,ipt3}(Ix, Ix)}), tt_tol);
            if isempty(Agttcur)
              Agttcur=curterm;
            else
              Agttcur = round(Agttcur + curterm, tt_tol);
            end
          end
        end
      end

      % Add the current TT matrix to the total
      if isempty(Agtt)
        Agtt = round(Agttcur, tt_tol);
      else
        Agtt = round(Agtt + Agttcur, tt_tol);
      end
    end
  end
  % add time operator
  AD = round(tkron(tt_matrix(MX1t(It,It)),Agtt) ...
    + tkron(tt_matrix(MX2t(It,It)),Agtt),tt_tol);

  %% Build global operator
  Att = round(AT + AD,tt_tol);

  %% Get the rhs term
  LLtt= amen_cross_zero(Ctt, @(x) cross_fun_nD(x,gfn),tt_tol,'verb',0);
  MMxtt = matrices_to_tt_matrix_fn(repmat({ML(Ix,:)}, 1, d));
  MMtt = tkron(tt_matrix(MLt(It,:)),MMxtt);
  gtt = amen_mv(MMtt,LLtt,tt_tol,'verb',0);

  %% Forming root finding problem

  MMxtt = matrices_to_tt_matrix_fn(repmat({ML(Ix,Ix)}, 1, d));
  MMtt2 = tkron(tt_matrix(MLt(It,It)),MMxtt);

  % F = @(y) round(amen_mv(Att,y,tt_tol,'verb',0) ...
  %   - amen_mv(MMtt2,y + y.^3,tt_tol,'verb',0) - gtt, tt_tol);
  % JF = @(y) round(Att - MMtt2 + ...
  %   3*amen_mm(MMtt2,make_tt_to_operator(y.^2),tt_tol),tt_tol);

  F = @(y) round(amen_mv(Att,y,tt_tol,'verb',0) ...
    - amen_mv(MMtt2, round(y - y.^3,tt_tol),tt_tol,'verb',0) - gtt, tt_tol);
  JF = @(y) round(Att - 2*round(MMtt2.*make_tt_to_operator(y - y.^3),tt_tol),tt_tol);

  %% Newton
  epsk = 1e-2;
  U0 = tt_zeros(size(gtt));
  Newton_eps = max(tt_tol,1e-6);
  %%
  nrmFU0 = norm(F(U0));
  maxiter = 100;
  du = U0;
  for k = 1:maxiter
    tsiter = datetime;
    % epsk = epsk*epsfactor;
    if 1
      fprintf('\n******** Newton iter = %d ************\n\n', k);
      fprintf('epsk = %.5e \n', epsk);
    end

    % compute Jacobian with JFbc
    Jmatk = JF(U0);

    % use GMRES to solve for du
    b = -F(U0);

    du = amen_solve2(Jmatk,b,epsk,'nswp',30,'verb',0);

    % tegmres = datetime;
    % fprintf('Linear Solve Time = %.2f\n', seconds(tegmres-tsiter));

    %update U1
    for alpha = [1,1/2,1/4,1/8,1/16]
      U1 = round(U0 + alpha*du, max(epsk,tt_tol));
      if norm(F(U1))<norm(F(U0))
        break;
      end
    end

    %compute local error
    local_err=norm(U1-U0)/norm(U0);
    res_norm = norm(F(U1))/nrmFU0;

    % #* set epsk
    epsk = min(min([local_err,res_norm]),epsk);
    % epsk = min([local_err,res_norm]);
    if 1
      fprintf('apha = %.2e \n', alpha);
      fprintf('u_err = %.5e,  Fu_ratio = %.5e \n', local_err, res_norm)
      fprintf('iter time = %.2f \n', seconds(datetime-tsiter));
    end

    if (local_err<Newton_eps) || (res_norm<Newton_eps)
      break;
    end
    U0 = U1;
  end
  utt = U1;

  %%
  tt_time = seconds(datetime-ts);
  Agcomp = compress_ratio_tt(Att);
  ucomp = compress_ratio_tt(utt);
  %% compute error
  uexacttt= amen_cross_zero(Ctt, @(x) cross_fun_nD(x,uexactfn),tt_tol,'verb',0);
  uexacttt = tt_get_inner(uexacttt,{2:Nt,2:Nx-1,2:Nx-1,2:Nx-1});

  Errtt(jjns)=norm(utt-uexacttt)/norm(uexacttt);
  utrunccomp = compress_ratio_tt(round(utt,Errtt(jjns)));
  
  %% store the results
  c.NewtonIter = k;
  c.hx = hx;
  c.ht = ht;
  c.tt_tol = tt_tol;
  c.error = Errtt(jjns);
  c.Agttcomp = Agcomp;
  c.Agttrank = Agtt.r;
  c.time = tt_time;

  R{jjns,1} = c;
  %% print out errors
  if jjns==1
    fprintf('Ex = %d, ',Ex);
    fprintf('qtt error = %.2e \n',Errtt(jjns));
  else
    fprintf('Ex = %d, tt Err = %.5e , convrate = %.5f\n',Ex,Errtt(jjns),...
      ( log(Errtt(jjns)) - log(Errtt(jjns-1)) )...
      /log((NEs(jjns-1)-1)/(NEs(jjns)-1)));
  end
  fprintf('hx = %.2e - tt tol = %.2e \n', hx, tt_tol);
  fprintf('Elapsed Time = %.5f seconds \n',tt_time)
  fprintf('Ag compress = %.2e \n', Agcomp);
  fprintf('u compress = %.2e \n', ucomp);
  fprintf('truncated u compress = %.2e \n', utrunccomp);

  %% save
  save(fname,'NEs','uexactfn','gfn','R');
  fprintf('Result is saved for Nx = %d  in file %s \n', Nx, fname);
end

