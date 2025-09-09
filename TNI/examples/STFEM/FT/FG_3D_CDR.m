close all; clear; clc;

addpath(genpath('../../../matlab/STFEM/src/'))
addpath(genpath('../../../matlab/utils/chebfun/'))
addpath(genpath('../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../matlab/utils/ttfunc/'))

%% problem definitions
uexactfn = @(t,x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z).*sin(pi*t) ...
  + sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z).*sin(2*pi*t)...
  + sin(3*pi*x).*sin(3*pi*y).*sin(3*pi*z).*sin(3*pi*t);

cfn = @(x,y,z) exp(-x -y -z);

afn = @(x,y,z) cos(pi*x).*cos(pi*y).*cos(pi*z) + 1;
bfn = {@(x,y,z) x, @(x,y,z) y, @(x,y,z) z};

% compute terms
dudt = compute_dudt_fn(uexactfn);
Dfn = compute_3D_diffusion_fn(uexactfn,afn);
Rfn = compute_3D_reaction_fn(uexactfn,cfn);
Cvecfn = compute_3D_convection_fn(uexactfn, bfn);

f = @(t,x,y,z) dudt(t,x,y,z) + Dfn(t,x,y,z) + Cvecfn(t,x,y,z) + Rfn(t,x,y,z);
%%
d = 3;
NEs = 1 + 2.^[2]; %full run [2:4]
%%
R = cell(numel(NEs),1);
fname='../plot_data/FG_3D_CDR.mat';

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

  % tt_tol = 1e-8;
  tt_tol = 0.1*hx^2;
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
  NNL = sparse(Ex, 2*Ex);  % Using sparse matrix for efficiency
  NNLt = sparse(Et, 2*Et);  % Using sparse matrix for efficiency

  % Populate NNL matrix
  for ii=2:Ex
    NNL(ii,ii+(ii-2))=1; NNL(ii,ii+(ii-2)+1)=1;
  end
  NNL(1,1)=1; NNL(Nx, 2*Ex)=1;

  NNR=NNL';

    % Populate NNL matrix
  for ii=2:Et
    NNLt(ii,ii+(ii-2))=1; NNLt(ii,ii+(ii-2)+1)=1;
  end
  NNLt(1,1)=1; NNLt(Nt, 2*Et)=1;

  NNRt=NNLt';


  % Assemble stiffness matrix and load vector
  [A1,A2,M1,M2,B1,B2] = nonlinear_Mat(X(1:2));
  MT = Mat_Time(T(1:2));

  AA = {kron(eye(Ex),A1),kron(eye(Ex),A2)};
  MM = {kron(eye(Ex),M1),kron(eye(Ex),M2)};
  BB = {kron(eye(Ex),B1),kron(eye(Ex),B2)};
  % Element stiffness matrix for a Q1 element in 2D
  Bs=(hx)*(1/6)*[2 1;1 2];
  ML = NNL*kron(eye(Ex),Bs)*NNR;
  MX1 = NNL*kron(eye(Ex),M1)*NNR;
  MX2 = NNL*kron(eye(Ex),M2)*NNR;
  

  % matrix for time dimension
  ZT = NNLt*kron(eye(Et),MT)*NNRt;
  MLt = NNLt*kron(eye(Et),Bs)*NNRt;

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

      %% Build 1D matrix operator
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
  AD = round(tkron(tt_matrix(MX1(It,It)),Agtt) ...
    + tkron(tt_matrix(MX2(It,It)),Agtt),tt_tol);

  %% Build Convection operator - x term
  btt = amen_cross_zero(Cxtt, @(x) cross_fun_nD(x,bfn{1}),tt_tol,'verb',0);
  G = core2cell(btt);
  nt1 = btt.r(2);
  nt2 = btt.r(3);
  
  Agtt = [];
  for j1 = 1:nt1
    for j2 = 1:nt2

      %% Build 1D matrix operator
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
          AB{idim,ipt} = NNL * (BB{ipt} * diag_CC{idim,ipt}) * NNR;
          MG{idim,ipt} = NNL * (MM{ipt} * diag_CC{idim,ipt}) * NNR;
        end
      end

      % Calculate the current TT matrix
      Agttcur =[];
      for ipt = 1:2
        for ipt2 = 1:2
          for ipt3 = 1:2
            curterm = round(...
              matrices_to_tt_matrix_fn({AB{1,ipt}(Ix, Ix), ...
              MG{2,ipt2}(Ix, Ix), MG{3,ipt3}(Ix, Ix)}), tt_tol);
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
  ACxtt = Agtt;

  %% Build Convection operator - y term
  btt = amen_cross_zero(Cxtt, @(x) cross_fun_nD(x,bfn{2}),tt_tol,'verb',0);
  G = core2cell(btt);
  nt1 = btt.r(2);
  nt2 = btt.r(3);
  
  Agtt = [];
  for j1 = 1:nt1
    for j2 = 1:nt2

      %% Build 1D matrix operator
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
          AB{idim,ipt} = NNL * (BB{ipt} * diag_CC{idim,ipt}) * NNR;
          MG{idim,ipt} = NNL * (MM{ipt} * diag_CC{idim,ipt}) * NNR;
        end
      end

      % Calculate the current TT matrix
      Agttcur =[];
      for ipt = 1:2
        for ipt2 = 1:2
          for ipt3 = 1:2
            curterm = round(...
              matrices_to_tt_matrix_fn({MG{1,ipt}(Ix, Ix), ...
              AB{2,ipt2}(Ix, Ix), MG{3,ipt3}(Ix, Ix)}), tt_tol);
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
  ACytt = Agtt;
  %% Build Convection operator - z term
  btt = amen_cross_zero(Cxtt, @(x) cross_fun_nD(x,bfn{3}),tt_tol,'verb',0);
  G = core2cell(btt);
  nt1 = btt.r(2);
  nt2 = btt.r(3);
  
  Agtt = [];
  for j1 = 1:nt1
    for j2 = 1:nt2

      %% Build 1D matrix operator
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
          AB{idim,ipt} = NNL * (BB{ipt} * diag_CC{idim,ipt}) * NNR;
          MG{idim,ipt} = NNL * (MM{ipt} * diag_CC{idim,ipt}) * NNR;
        end
      end

      % Calculate the current TT matrix
      Agttcur =[];
      for ipt = 1:2
        for ipt2 = 1:2
          for ipt3 = 1:2
            curterm = round(...
              matrices_to_tt_matrix_fn({MG{1,ipt}(Ix, Ix), ...
              MG{2,ipt2}(Ix, Ix), AB{3,ipt3}(Ix, Ix)}), tt_tol);
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
  ACztt = Agtt;
  %% 
  ACspace = ACxtt + ACytt + ACztt;
  % add time operator
  AC = round(tkron(tt_matrix(MX1(It,It)),ACspace) ...
    + tkron(tt_matrix(MX2(It,It)),ACspace),tt_tol);
  
  %% Reaction Operator
  ctt= amen_cross_zero(Cxtt, @(x) cross_fun_nD(x,cfn),tt_tol,'verb',0);
  G = core2cell(ctt);
  nt1 = ctt.r(2);
  nt2 = ctt.r(3);
  
  Agtt = [];
  for j1 = 1:nt1
    for j2 = 1:nt2

      %% Build 1D matrix operator
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
          MG{idim,ipt} = NNL * (MM{ipt} * diag_CC{idim,ipt}) * NNR;
        end
      end

      % Calculate the current TT matrix
      Agttcur =[];
      for ipt = 1:2
        for ipt2 = 1:2
          for ipt3 = 1:2
            curterm = round(...
              matrices_to_tt_matrix_fn({MG{1,ipt}(Ix, Ix), ...
              MG{2,ipt2}(Ix, Ix), MG{3,ipt3}(Ix, Ix)}), tt_tol);
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
  AR = round(tkron(tt_matrix(MX1(It,It)),Agtt) ...
    + tkron(tt_matrix(MX2(It,It)),Agtt),tt_tol);
  %% %%%%%%%%%%%%% Build global %%%%%%%%%%%%%%% %%
  Afg = full(AT) + full(AC) + full(AD) + full(AR);
  
  %% Get the rhs term
  LLtt= amen_cross_zero(Ctt, @(x) cross_fun_nD(x,f),tt_tol,'verb',0);
  MMxtt = matrices_to_tt_matrix_fn(repmat({ML(Ix,:)}, 1, d));
  MMtt = tkron(tt_matrix(MLt(It,:)),MMxtt);
  F_newtt = amen_mv(MMtt,LLtt,tt_tol,'verb',0);
  F_new = full(F_newtt);

  %% linear solve
  u = Afg\F_new;
  
  tt_time = seconds(datetime-ts);
  
  %% compute error
  uexacttt= amen_cross_zero(Ctt, @(x) cross_fun_nD(x,uexactfn),tt_tol,'verb',0);
  uexacttt = tt_get_inner(uexacttt,{2:Nt,2:Nx-1,2:Nx-1,2:Nx-1});
  uexact = full(uexacttt);
  
  Errtt(jjns)=norm(u-uexact)/norm(uexact);

  %% store the results
  c.hx = hx;
  c.ht = ht;
  c.error = Errtt(jjns);
  c.time = tt_time;
  R{jjns} = c;
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

  %% save
  save(fname,'NEs','uexactfn','afn','bfn','cfn','f','R');
  fprintf('Result is saved for Nx = %d \n', Nx);

end

