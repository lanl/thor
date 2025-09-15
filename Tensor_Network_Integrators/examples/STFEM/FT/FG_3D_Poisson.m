close all; clear; clc;

addpath(genpath('../../../matlab/STFEM/src/'))
addpath(genpath('../../../matlab/utils/chebfun/'))
addpath(genpath('../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../matlab/utils/ttfunc/'))

%% problem definitions
uexactfn = @(x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z);
afn = @(x,y,z) 1 + cos(pi*(x+y)).*cos(pi*z);
f = compute_3D_rhs_fn(uexactfn,afn);
d = 3;
NEs = 1 + 2.^[3:4]; %full run 3:5

amen_cross = @amen_cross_zero;
%%
R = cell(numel(NEs),1);
fname='../plot_data/FG_3D_Poisson.mat';

for jjns = [1,1,1:numel(NEs)]
  %% parameters
  tic;

  Lx = 1; % Length of the domain in x-direction
  Ly = 1; % Length of the domain in y-direction
  Lz = 1;
  Ex = NEs(jjns); % Number of elements in x-direction
  hx = Lx / Ex; % Element size in x-direction
  Nx = Ex + 1; % Number of nodes in x-direction

  tt_tol = 1e-6;

  %% Create mesh
  X=0:hx:1;

  %% create grid in tt
  Itt = repmat({ones(1, Nx)}, 1, d);
  Ctt = cell(1, d);
  for ic = 1:d
    temp = Itt;
    temp{ic} = X;
    Ctt{ic} = cell2core(tt_tensor,temp);
  end

  %% %%%%%%%%%%% Construction of the left hand side matrix
  % Initialize NNL with zeros
  NNL = sparse(Ex, 2*Ex);  % Using sparse matrix for efficiency
  % Populate NNL matrix
  for ii=2:Ex
    NNL(ii,ii+(ii-2))=1; NNL(ii,ii+(ii-2)+1)=1;
  end
  NNL(1,1)=1; NNL(Nx, 2*Ex)=1;
  NNR=NNL';

  % Assemble stiffness matrix and load vector
  [A1,A2,M1,M2] = nonlinear_Mat(X(1:2));
  AA{1} = kron(eye(Ex),A1);
  AA{2} = kron(eye(Ex),A2);
  MM{1} = kron(eye(Ex),M1);
  MM{2} = kron(eye(Ex),M2);

  % Element stiffness matrix for a Q1 element in 2D
  Bs=(hx)*(1/6)*[2 1;1 2];
  ML = NNL*kron(eye(Ex),Bs)*NNR;
  %%  FG Tensor Global Assemble
  % nterm = 2; %number of rank-1 term in the decomposition of funciton a
  att= amen_cross(Ctt, @(x) cross_fun_nD(x,afn),tt_tol,'verb',0);
  %%
  G = core2cell(att);
  nt1 = att.r(2);
  nt2 = att.r(3);

  c.build1Dtime = toc;
  %% Build the global operator for only interior nodes
  tic;
  Agtt = 0;
  Ix = [2:Ex]; % interior system

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
      Agttcur =0;
      for ipt = 1:2
        for ipt2 = 1:2
          for ipt3 = 1:2
            curterm = matrices_to_tt_matrix_fn({AG{1,ipt}(Ix, Ix), ...
              MG{2,ipt2}(Ix, Ix), MG{3,ipt3}(Ix, Ix)}) ...
              + matrices_to_tt_matrix_fn({MG{1,ipt}(Ix, Ix), ...
              AG{2,ipt2}(Ix, Ix),MG{3,ipt3}(Ix, Ix)}) ...
              + matrices_to_tt_matrix_fn({MG{1,ipt}(Ix, Ix), ...
              MG{2,ipt2}(Ix, Ix),AG{3,ipt3}(Ix, Ix)});
            
            Agttcur = Agttcur + full(curterm);

          end
        end
      end

      Agtt = Agtt + Agttcur;
    end
  end
  c.buildTTopstime = toc;
  %% Get the rhs term
  tic;
  LLtt= amen_cross_zero(Ctt, @(x) cross_fun_nD(x,f),tt_tol,'verb',0);
  LLtt = full(LLtt);
  MMtt = full(matrices_to_tt_matrix_fn(repmat({ML(Ix,:)}, 1, d)));
  % F_newtt = amen_mv(MMtt,LLtt,tt_tol,'verb',0);
  F_new = MMtt*LLtt;

  c.rhsbuildtime = toc;

  %% linear solve
  tic;
  uqtt = Agtt\F_new;

  tqTTsolve = toc;
  %% compute error
  uexacttt= amen_cross_zero(Ctt, @(x) cross_fun_nD(x,uexactfn),tt_tol,'verb',0);
  uexacttt = tt_get_inner(uexacttt,repmat({2:Nx-1},1,d));
  uexacttt = tt_reshape(uexacttt, size(uqtt),tt_tol);

  Errtt(jjns)=norm(uqtt-full(uexacttt))/norm(uexacttt);

  %% store the results
  c.hx = hx;
  c.tt_tol = tt_tol;
  c.error = Errtt(jjns);
  c.buildtime = c.build1Dtime + c.buildTTopstime + c.rhsbuildtime;
  c.TTsolvetime = tqTTsolve;
  c.time = c.buildtime+c.TTsolvetime;

  R{jjns,1} = c;
  %% print out errors
  fprintf('%s\n',repmat('*',30,1));
  if jjns==1
    fprintf('Ex = %d, ',Ex);
    fprintf('qtt error = %.2e \n',Errtt(jjns));
  else
    fprintf('Ex = %d, qtt Err = %.5e , convrate = %.5f\n',Ex,Errtt(jjns),...
      ( log(Errtt(jjns)) - log(Errtt(jjns-1)) )...
      /log((NEs(jjns-1)-1)/(NEs(jjns)-1)));
  end
  fprintf('hx = %.2e - tt tol = %.2e \n', hx, tt_tol);
  fprintf('Elapsed time = %.2f s \n',c.time);

  %% save
  save(fname,'NEs','uexactfn','afn','f','R');
  fprintf('Result is saved for Nx = %d \n', Nx);
end

