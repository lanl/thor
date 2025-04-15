close all; clear; clc;
run setup1DSlab.m;
%% rescale the problem
tic;
a = 0;
b = 1;
alpha = 0; % no second term
% scatterting parameters
sigma_s = importdata('./Data/sigma_s_8groups.txt');
% fission parameters
nusigma_f = load('./Data/eig_nusigma_f.txt');
chi = load('./Data/eig_chi.txt');

% Energy parameters
% sigma = 0.45391;
sigma =[ 2.3336e-01,2.3105e-01, 3.0084e-01, 2.8809e-01, ...
  3.8538e-01, 5.0120e-01, 6.1661e-01, 7.8866e-01 ];

nE = numel(sigma); % index g
nmu = 128; % index l
nx = 128; % index i - number of EDGES
dx = (b-a)/(nx-1); %spacial step

% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]);
mu = flip(mu');
wgt = flip(wgt'/2);

%% forming the system in the tensor format
% positive mu part
Lposmumat = 1/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));
Lposgmat = diag(ones(nx,1)) + diag(ones(nx-1,1),-1);

% negative mu part
Lnegmumat = 1/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
Lneggmat = diag(ones(nx,1)) + diag(ones(nx-1,1),1);

%% H1

G = cell(1,3);
G{1} = reshape(Lposmumat,[1,size(Lposmumat),1]);
mupos = diag([mu(1:nmu/2), zeros(1,nmu/2)]);
G{2} = reshape(mupos,[1,nmu,nmu,1]);
G{3} = reshape(eye(nE),[1,size(eye(nE)),1]);
ttLposmu = cell2core(tt_matrix, G);

% Lposg
G{1} = reshape(Lposgmat,[1,size(Lposgmat),1]);
G{2} = G{2}~=0;
G{3} = reshape(diag(sigma/2),[1,nE,nE,1]);
ttLposg = cell2core(tt_matrix, G);
ttH1 = ttLposmu + ttLposg;

%% H2

G = cell(1,3);
G{1} = reshape(Lnegmumat,[1,size(Lnegmumat),1]);
muneg = diag([zeros(1,nmu/2), mu(nmu/2+1:nmu)]);
G{2} = reshape(muneg,[1,nmu,nmu,1]);
G{3} = reshape(eye(nE),[1,size(eye(nE)),1]);
ttLnegmu = cell2core(tt_matrix, G);

% only update G{1} and G{3}
G{1} = reshape(Lneggmat,[1,size(Lneggmat),1]);
G{2} = G{2}~=0;
G{3} = reshape(diag(sigma/2),[1,nE,nE,1]);
ttLnegg = cell2core(tt_matrix, G);

ttH2 = ttLnegmu + ttLnegg;

%forming H in tt
ttH = ttH1 + ttH2;


%% Construct the scattering operator Hs

% Hs positive
tempXpos = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1));
tempXpos_noBC = tempXpos;
tempXpos_noBC(1,:) = 0;

tempmupos = [ones(nmu/2,1);zeros(nmu/2,1)]*wgt';

G = cell(1,3);
G{1} = reshape(tempXpos_noBC,[1,nx,nx,1]);
G{2} = reshape(tempmupos,[1,nmu,nmu,1]);
G{3} = reshape(sigma_s,[1,nE,nE,1]);
ttHspos = cell2core(tt_matrix,G);

% Hs negative
tempXneg = permute(tempXpos,[2,1]);
tempXneg_noBC = tempXneg;
tempXneg_noBC(end,:) = 0;

G = cell(1,3);
G{1} = reshape(tempXneg_noBC,[1,nx,nx,1]);
G{2} = reshape(flipud(tempmupos),[1,nmu,nmu,1]);
G{3} = reshape(sigma_s,[1,nE,nE,1]);
ttHsneg = cell2core(tt_matrix,G);
ttHs = ttHspos + ttHsneg;


%% Construct the fission Hf
G = core2cell(ttHspos);
G{3} = reshape(chi.*nusigma_f',[1,nE,nE,1]);
ttHfpos = cell2core(tt_matrix, G);

G = core2cell(ttHsneg);
G{3} = reshape(chi.*nusigma_f',[1,nE,nE,1]);
ttHfneg = cell2core(tt_matrix, G);

ttHf = ttHfpos + ttHfneg;


%% Fixed-point scheme in tensor train
tol = 1e-6;
% fixed point schemes
ttPsi0 = tt_rand_pos([nx;nmu;nE],3,1);
k = 1;
niter = 100;
convrate = nan(1,niter);

for iter = 1:niter
  fprintf('Iteration = %d/%d \n convrate = %.5e \n', ...
    iter,niter, convrate(max(iter-1,1)));
  RHS = amen_mv(ttHs,ttPsi0,tol) + 1/k*amen_mv(ttHf,ttPsi0,tol);

  %  Solve a linear system for Psi1 = Htt\RHS;
  ttPsi1 = amen_solve2(ttH,RHS,tol);

  %update eigenvalue
  k = k*sum(amen_mv(ttHf,ttPsi1,tol))./sum(amen_mv(ttHf,ttPsi0,tol));
  %   fprintf('iter = %d - ||Psi1 - Psi0|| = %.5e \n',iter, norm(ttPsi1-ttPsi0));
  convrate(iter) = norm(ttPsi1-ttPsi0)/norm(ttPsi0);
  if convrate(iter) < tol
    fprintf('The fixed point scheme converged at iter = %d \n', iter)
    break
  end
  ttPsi0 = ttPsi1;
end
%normalize ttPsi1
ttPsi1 = ttPsi1/norm(ttPsi1);
% %%

fprintf('-------> k = %.5e \n',k);
fprintf('Psi1 compress ratio = %.5e \n ', compress_ratio_tt(ttPsi1));
toc
