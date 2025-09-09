close all; clear; clc;
run setup1DSlab.m;

%% problem setup
% solution Psi (261 Energy - 256 mu - 256 x)

%% load parameters
Param.sigma_t = importdata('./Data/sigma_t.txt');
Param.nusigma_f = importdata('./Data/nusigma_f.txt');
Param.vel = importdata('./Data/vel.txt');
Param.sigma_s = importdata('./Data/sigma_s.txt');
Param.mu = importdata('./Data/mu.txt');
Param.wgt = importdata('./Data/wgt.txt');
Param.chi = importdata('./Data/chi.txt');

convtol = 1e-6;
ConvRate =[];
nIter = 100; % fixed point iteration
%% rescale the problem
alpha = 0; % no second term
nE = 261; % index g
nmu = 8; % index l
nx = 9; % index i
dx = 4.2/(nx-1); %space step

% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]);
mu = flip(mu);
wgt = flip(wgt'/2);

% Energy parameters
sigma_t = Param.sigma_t;
nusigma_f = Param.nusigma_f;
sigma_s = Param.sigma_s;
chi = Param.chi;


%% %%%%%%%%% Fix point scheme  %%%%%%%%%%%%%%%
%{


%}
%% Initial guess - random rank 1 tensor
% tPsi0 = tt_rand([nE,nmu,nx],3,1);
% Psi0 = tensor(full(tPsi0, size(tPsi0))); %dim nE x nmu x nx
Psi0 = tensor(rand(nE,nmu,nx));
% set zero boundary in the x direction
Psi0(:,:,1) = 0;
Psi0(:,:,nx) = 0;

for iter = 1:nIter
  %% compute RHS0
  % Hmu
  mu_mat = diag(mu)/dx;

  Hmu = ttt(tensor(eye(nE)),ttt(tensor(mu_mat),  tensor(eye(nx-2))));
  
  Hg = ttt(tensor(1/2*diag(sigma_t)), ttt(tensor(eye(nmu)),tensor(eye(nx-2))));

  RHS0 = ttt((Hmu-Hg),Psi0(:,:,2:end-1),[2,4,6],[1,2,3]);
  %% compute the Lterm
  %tempL can be used for both RHS1 and RHS2
  Lterm = 0.5*(Psi0(:,:,2:end-1) + Psi0(:,:,1:end-2)); %tempL dim nE x nmu x (nx-2)

  %contract along mu dimension
  Lterm = tensor(Lterm);
  %ttv is tensor time vector
  Lterm = ttv(Lterm,wgt,2); % dimension nE x (nx-2)


  %% compute the RHS1
  RHS1 = 0.5*ttv(Lterm,nusigma_f,1); % contract along the Energy dimension
  % form RHS1 as 3D tensor
  % outer product with ones in mu direction
  RHS1 = ttt(tensor(ones(1,nmu),nmu), RHS1);
  % outer product with chi/2
  RHS1 = ttt(tensor(chi/2,nE), RHS1); % dimension of RHS1 nE x nmu x (nx-2);

  %% compute the RHS2
  RHS2 = 0.5*ttm(Lterm,sigma_s,1); % multiply sigma_s(g,g') with Lterm(g,i)
  RHS2 = ttt(tensor(ones(1,nmu),nmu), RHS2);
  RHS2 = permute(RHS2,[2,1,3]);

  %% solve a linear system for Psi1
  tempA = reshape(permute((Hmu+Hg),[1,3,5,2,4,6]),[nE*nmu*(nx-2),nE*nmu*(nx-2)]);
  tempb = RHS0 + RHS1 + RHS2;
  tempb = tempb(:);
  tempx = mldivide(double(tempA),double(tempb)); %solve the linear system
  % reshape to new Psi
  Psi1 = reshape(tempx,[nE, nmu, nx-2]);

  %% boundary treatmeant
  Psi1 = cat(3,zeros(nE,nmu,1),Psi1,zeros(nE,nmu,1));%concatenate zero layers
  Psi1 = tensor(Psi1);
  %% convergence error
  converr = norm(Psi1-Psi0)/norm(Psi0);
  ConvRate = [ConvRate; converr];
  fprintf('iter = %3.d, converr = %.5e \n', iter, converr);
  if converr<=convtol
    break;
  else
    Psi0 = Psi1;
  end

end
