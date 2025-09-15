close all; clear; clc;
run setup1DSlab.m;

%% problem setup
%{ 
solution Psi (2 mu - 3 x)
mu1*dPsi/dx + sigma*Psi = Q1
mu2*dPsi/dx + sigma*Psi = Q2
%}

data.sigma_s = importdata('./Data/sigma_s_8groups.txt');

%% rescale the problem
a = 0;
b = 1;
alpha = 0; % no second term
% Energy parameters

% sigma = 0.45391;

sigma =[ 2.3336e-01,2.3105e-01, 3.0084e-01, 2.8809e-01, ...
  3.8538e-01, 5.0120e-01, 6.1661e-01, 7.8866e-01 ];

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
Lposmu = ttt(tensor(Lposmumat),tensor(eye(nE)));
Lposgmat = diag(ones(nx,1)) + diag(ones(nx-1,1),-1);
Lposg = ttt(tensor(Lposgmat), tensor(diag(sigma/2)));
Lpos = Lposmu+Lposg;

% negative mu part
Lnegmumat = mu(2)/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
Lnegmu = ttt(tensor(Lnegmumat), tensor(eye(nE)));

Lneggmat = diag(ones(nx,1)) + diag(ones(nx-1,1),1);
Lnegg = ttt(tensor(Lneggmat), tensor(diag(sigma/2)));
Lneg = Lnegmu + Lnegg;

%%
Lpos = tensor(Lpos);
Lneg = tensor(Lneg);
H1 = ttt(tensor(Lpos), tensor([1,0;0,0]));
H1 = permute(H1,[1,2,5,6,3,4]);
H2 = ttt(tensor(Lneg), tensor([0,0;0,1]));
H2 = permute(H2,[1,2,5,6,3,4]);
H = H1 + H2;

%% %%% tensorize H1
if 1
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
  %check error
  H1mat = convert_ten_to_mat_fn(H1,nx,nmu,nE);
  fprintf('ttH1 error = %.5e \n', norm(full(ttH1)-H1mat));
end

%% %%% tensorize H2 
if 1
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
  %check error
  H2mat = convert_ten_to_mat_fn(H2,nx,nmu,nE);
  fprintf('ttH2 error = %.5e \n', norm(full(ttH2)-H2mat));
  
  %forming H in tt
  ttH = ttH1 + ttH2;
end

%% Construct the scattering operator tensor
% it has 2 parts for neg/pos mu
% no need to include the boundary condition here
sigma_s = data.sigma_s;

tempXpos = tensor(1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1)));
tempXpos_noBC = tempXpos;
tempXpos_noBC(1,:) = 0;
tempmupos = tensor([1;0]*wgt');
Hspos = ttt(tempXpos_noBC,tempmupos);
% Hsposmat = reshape(permute(Hspos,[1,3,2,4]),[nx*nmu,nx*nmu]);

tempXneg = permute(tempXpos,[2,1]);
tempXneg_noBC = tempXneg;
tempXneg_noBC(end,:) = 0;

Hsneg = ttt(tempXneg_noBC, flipud(tempmupos));
% Hsnegmat = reshape(permute(Hsneg,[1,3,2,4]),[nx*nmu,nx*nmu]);

%
Hs = ttt(Hspos + Hsneg, tensor(sigma_s));
%% %%% tensorize Hs %%%%%%
if 1
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
%   Hsposmat = convert_ten_to_mat_fn(Hspos,nx, nmu, nE);
%   fprintf('tt - Hspos error = %.5e \n', norm(full(ttHspos)-Hsposmat));

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
  Hsmat = convert_ten_to_mat_fn(Hs,nx, nmu, nE);
  fprintf('tt - Hs error = %.5e \n', norm(full(ttHs)-Hsmat));
end

%% Construct the fission operator Hf
% resuse Hspos and Hsneg, which is the interpolation and mu_quadrature term
nusigma_f = load('./Data/eig_nusigma_f.txt');
chi = load('./Data/eig_chi.txt');
Hf = ttt(Hspos + Hsneg, tensor(chi.*nusigma_f'));

%% %%%%%%% tensorize Hf %%%%%%%%
G = core2cell(ttHspos);
G{3} = reshape(chi.*nusigma_f',[1,nE,nE,1]);
ttHfpos = cell2core(tt_matrix, G);

G = core2cell(ttHsneg);
G{3} = reshape(chi.*nusigma_f',[1,nE,nE,1]);
ttHfneg = cell2core(tt_matrix, G);

ttHf = ttHfpos + ttHfneg;

Hfmat = convert_ten_to_mat_fn(Hf, nx, nmu, nE);
fprintf('tt - Hf error = %.5e \n', norm(full(ttHf)-Hfmat));


%% accumulate to H
A = H - Hs;
B = Hf;
%% solve the eigenvalue problem
% Amat = reshape(permute(A,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]);
% Bmat = reshape(permute(B,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]);
% % solve the eigvenvalue problem
% [V,D] = eig(double(Amat), double(Bmat));

% [V,D] = eig(inv(double(Bmat) + 1e-6*eye(size(Bmat)))*double(Amat));
%% Fixed-point scheme in fullgrid
Hsmat = convert_ten_to_mat_fn(Hs,nx,nmu,nE);
Hmat = convert_ten_to_mat_fn(H,nx,nmu,nE);
Hfmat = convert_ten_to_mat_fn(Hf,nx,nmu,nE);
% Psi0 = rand(nx*nmu*nE,1);
Psi0tt = tt_rand_pos([nx;nmu;nE],3,1);
Psi0 = full(Psi0tt);
k = 1;
for iter = 1:100
  RHS = Hsmat*Psi0 + 1/k*Hfmat*Psi0;

  Psi1 = Hmat\RHS;
  k = k*sum(Hfmat*Psi1)/sum(Hfmat*Psi0);
  fprintf('||Psi1-Psi0|| = .%5e \n',norm(Psi1-Psi0));
  Psi0 = Psi1;
end
Psi1 = Psi1/norm(Psi1);
fprintf('fixed point k = %.5e \n', k);
%% solve by tensor train
% H_for_tt = permute(double(H),[1,3,5,2,4,6]);
% Htt = tt_matrix(H_for_tt);
Htt = ttH;
Hstt = ttHs;
Hftt = ttHf;

tol = 1e-6;
% fixed point schemes
Psi0tt = tt_rand_pos([nx;nmu;nE],3,1);
ktt = 1;
for iter = 1:100
  RHS = amen_mv(Hstt,Psi0tt,tol) + 1/ktt*amen_mv(Hftt,Psi0tt,tol);

  %  Solve a linear system for Psi1 = Htt\RHS;
  Psi1tt = amen_solve2(Htt,RHS,tol);
  
  %update eigenvalue
  ktt = ktt*sum(amen_mv(Hftt,Psi1tt,tol))./sum(amen_mv(Hftt,Psi0tt,tol));
  fprintf('iter = %d - ||Psi1tt-Psi0tt|| = %.5e \n',iter, norm(Psi1tt-Psi0tt));
  Psi0tt = Psi1tt;
end
Psi1tt = Psi1tt/norm(Psi1tt);
fprintf('fixed point in TT - ktt = %.5e \n', ktt);
%%
fprintf('||k-ktt|| = %.5e \n', abs(k-ktt));
fprintf('||Psi1tt-Psi1||/||Psi1|| = %.5e \n', check_tt_rel_error(Psi1,Psi1tt));
fprintf('Psi1 compress ratio = %.5e \n ', compress_ratio_tt(Psi1tt));
