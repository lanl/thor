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
mu = flip(mu');
wgt = flip(wgt'/2);

%% forming the system in the tensor format

Lposmu = mu(1)/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));

Lposmu = ttt(tensor(Lposmu),tensor(eye(nE)));

Lposg = diag(ones(nx,1)) + diag(ones(nx-1,1),-1);
Lposg = ttt(tensor(Lposg), tensor(diag(sigma/2)));
Lpos = Lposmu+Lposg;

Lnegmu = mu(2)/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
Lnegmu = ttt(tensor(Lnegmu), tensor(eye(nE)));
Lnegg = diag(ones(nx,1)) + diag(ones(nx-1,1),1);
Lnegg = ttt(tensor(Lnegg), tensor(diag(sigma/2)));
Lneg = Lnegmu + Lnegg;

Lpos = tensor(Lpos);
Lneg = tensor(Lneg);

%%
H1 = ttt(tensor(Lpos), tensor([1,0;0,0]));
H2 = ttt(tensor(Lneg), tensor([0,0;0,1]));
H = H1+H2;
H = permute(H,[1,2,5,6,3,4]);
% current shape of H nx x nx x nE x nE x nmu x nmu

% adjust the boundary with coefficients c and d
% Q(1,1,:) = Q(1,1)/c(1);
% Q(nx,nmu,:) = Q(nx,nmu)/d(nmu);

%% Construct the scattering operator tensor
% it has 2 parts for neg/pos mu
% no need to include the boundary condition here
sigma_s = data.sigma_s;
% Hfp = ttt(tensor(diag([1,ones(1,nx-1)])),tensor([wgt(2)/2,0;0,0]));
% Hfp = ttt(Hfp,tensor(sigma_s/2));
% Hfm = ttt(tensor(diag([ones(1,nx-1),1])),tensor([0,0;0,wgt(1)/2]));
% Hfm = ttt(Hfm,tensor(sigma_s/2));
% Hf = Hfp+Hfm;

tempXpos = tensor(1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1)));
tempXpos_noBC = tempXpos;
tempXpos_noBC(1,:) = 0;
tempmupos = tensor([1;0]*wgt');

Hspos = ttt(tempXpos_noBC,tempmupos);
% Hfneg = ttt(ttt(permute(tempXpos,[2,1]),flipud(tempmupos)),tensor(sigma_s(1)/2));
% Hfpos = ttt(tempXpos,tensor([wgt(1),0; 0,0])) + ...
%   ttt(tempXpos_noBC,tensor([0, wgt(2); 0,0]));

Hfposmat = reshape(permute(Hspos,[1,3,2,4]),[nx*nmu,nx*nmu]);

tempXneg = permute(tempXpos,[2,1]);
tempXneg_noBC = tempXneg;
tempXneg_noBC(end,:) = 0;

Hsneg = ttt(tempXneg_noBC, flipud(tempmupos));
% Hfneg = ttt(tempXneg_noBC, tensor([0,0;wgt(1),0]))+ ...
%   ttt(tempXneg, tensor([0,0;0,wgt(2)]));

Hfnegmat = reshape(permute(Hsneg,[1,3,2,4]),[nx*nmu,nx*nmu]);

%
Hs = ttt(Hspos + Hsneg, tensor(sigma_s));

%% Construct the fission operator Hf
% resuse Hspos and Hsneg, which is the interpolation and mu_quadrature term
nusigma_f = load('./Data/eig_nusigma_f.txt');
chi = load('./Data/eig_chi.txt');
Hf = ttt(Hspos + Hsneg, tensor(chi.*nusigma_f'));
Hfmat = reshape(permute(Hf,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]);
Hfref = importdata('./Data/SigmaS_matrix_8groups.txt');

%% accumulate to H
A = H - Hs;
B = Hf;
%% solve the eigenvalue problem
Amat = reshape(permute(A,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]);
Bmat = reshape(permute(B,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]);
% solve the eigvenvalue problem
[V,D] = eig(double(Amat), double(Bmat));

% [V,D] = eig(inv(double(Bmat) + 1e-6*eye(size(Bmat)))*double(Amat));
%% Fixed-point scheme in fullgrid
Hsmat = double(reshape(permute(Hs,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]));
Hmat = double(reshape(permute(H,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]));
Hfmat = double(reshape(permute(Hf,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]));
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
Htt = convert_tensor_operator_to_tt_matrix_fn(H);
Hstt = convert_tensor_operator_to_tt_matrix_fn(Hs);
Hftt = convert_tensor_operator_to_tt_matrix_fn(Hf);

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

%%

% fprintf('Error between Htt and H = %.5e \n', temperr)
% fprintf('compress ratio in Hs = %.2e \n', compress_ratio_tt(tt_tensor(double(Hs))));
% fprintf('compress ratio in H-Hs = %.2e \n',compress_ratio_tt(Htt,Htt.r, prod(size(Htt),2)));
% Qtt = tt_tensor(Q);
% 
% Psitt = amen_solve2(Htt, Qtt, 1e-6);
% solerr = norm(Psi(:) - full(Psitt));
% 
% fprintf('Difference between Psi and Psitt = %.5e \n', solerr);

function Htt = convert_tensor_operator_to_tt_matrix_fn(H)
H_for_tt = permute(double(H),[1,3,5,2,4,6]);
Htt = tt_matrix(H_for_tt);
end