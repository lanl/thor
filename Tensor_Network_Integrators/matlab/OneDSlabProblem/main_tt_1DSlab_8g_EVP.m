% close all; 
clear; clc;
run setup1DSlab.m;
%% rescale the problem
a = 0;
b = 1;
alpha = 0; % no second term
% scatterting parameters
sigma_s = importdata('./Data/sigma_s_8groups.txt');
% fission parameters
nusigma_f = load('./Data/eig_nusigma_f.txt');
chi = load('./Data/eig_chi.txt');

% Energy parameters
sigma = 0.45391;
%sigma =[ 2.3336e-01,2.3105e-01, 3.0084e-01, 2.8809e-01, ...
  % 3.8538e-01, 5.0120e-01, 6.1661e-01, 7.8866e-01 ];

nE = numel(sigma); % index g
nmu = 2; % index l
nx = 3; % index i - number of EDGES
dx = (b-a)/(nx-1); %spacial step

% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]);
mu = flip(mu);
wgt = flip(wgt'/2);

%% forming the system in the tensor format
% positive mu part
Lposmumat = mu(1)/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));
Lposgmat = diag(ones(nx,1)) + diag(ones(nx-1,1),-1);

% negative mu part
Lnegmumat = mu(2)/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
Lneggmat = diag(ones(nx,1)) + diag(ones(nx-1,1),1);

%% H1

G = cell(1,3);
G{1} = reshape(Lposmumat,[1,size(Lposmumat),1]);
G{2} = reshape([1,0;0,0],[1,2,2,1]);
G{3} = reshape(eye(nE),[1,size(eye(nE)),1]);
ttLposmu = cell2core(tt_matrix, G);

% only update G{1} and G{3}
G{1} = reshape(Lposgmat,[1,size(Lposgmat),1]);
G{3} = reshape(diag(sigma/2),[1,nE,nE,1]);
ttLposg = cell2core(tt_matrix, G);
ttH1 = ttLposmu + ttLposg;

%% H2

G = cell(1,3);
G{1} = reshape(Lnegmumat,[1,size(Lnegmumat),1]);
G{2} = reshape([0,0;0,1],[1,2,2,1]);
G{3} = reshape(eye(nE),[1,size(eye(nE)),1]);
ttLnegmu = cell2core(tt_matrix, G);

% only update G{1} and G{3}
G{1} = reshape(Lneggmat,[1,size(Lneggmat),1]);
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

tempmupos = [1;0]*wgt';

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
for iter = 1:100
  RHS = amen_mv(ttHs,ttPsi0,tol) + 1/k*amen_mv(ttHf,ttPsi0,tol);

  %  Solve a linear system for Psi1 = Htt\RHS;
  ttPsi1 = amen_solve2(ttH,RHS,tol);
  
  %update eigenvalue
  k = k*sum(amen_mv(ttHf,ttPsi1,tol))./sum(amen_mv(ttHf,ttPsi0,tol));
  
  %better update here use Rayleigh's quotient

  fprintf('iter = %d - ||Psi1 - Psi0|| = %.5e \n',iter, norm(ttPsi1-ttPsi0));
  ttPsi0 = ttPsi1;
end
ttPsi1 = ttPsi1/norm(ttPsi1);
fprintf('fixed point in TT - ktt = %.5e \n', k);
% %%
% fprintf('||k-ktt|| = %.5e \n', abs(k-k));
% fprintf('||Psi1tt-Psi1||/||Psi1|| = %.5e \n', check_tt_rel_error(Psi1,ttPsi1));
% fprintf('Psi1 compress ratio = %.5e \n ', compress_ratio_tt(ttPsi1));

%% export data
if 1
  GH = core2cell(ttH);
  GHs = core2cell(ttHs);
  GHf = core2cell(ttHf);
  true_eig_val = k;
  filename ='./Results/Data_for_1DSlab_261g_EVP.mat';
  save(filename,'GH','GHs','GHf','tol','nx','nmu','nE','true_eig_val');
  
end