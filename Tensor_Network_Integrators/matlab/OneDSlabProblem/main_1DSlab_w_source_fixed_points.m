close all; clear; clc;
run setup1DSlab.m;

%% problem setup
%{
solution Psi (2 mu - 3 x)
mu1*dPsi/dx + sigma*Psi = Q1
mu2*dPsi/dx + sigma*Psi = Q2
%}

%% rescale the problem
a = 0;
b = 1;
alpha = 0; % no second term
nE = 1; % index g
nmu = 2; % index l
nx = 3; % index i - number of EDGES
dx = (b-a)/(nx-1); %spacial step

% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]);
mu = flip(mu);
wgt = flip(wgt'/2);

% Energy parameters
sigma = 0.45391;
true_sol = [0, 1.209762;7.237695e-1, 7.237695e-1; 1.209762, 0];

% RHS1 Q
Q = 1*ones(nx,nmu); %this is a vector of ones
Q(1,1) = 0; %BC mu pos
Q(nx,nmu) = 0; % BC mu neg

%% forming the system in the tensor format
mu_mat = diag(mu)/dx;
Hmu = ttt(tensor(eye(nx)),tensor(mu_mat));
Hg = sigma/2.*ttt(tensor(eye(nx)), tensor(eye(nmu)));
H1 = Hmu + Hg;

%% Inital Guess
Psi0 = tensor(rand(nx,nmu));
% Psi0 = true_sol;
Psi0(1,1) = 0;
Psi0(nx,nmu) = 0;

convtol = 1e-6;
maxiter = 100;

for iter = 1%:maxiter
  % compute the RHS (shift term, Q terms and boundary terms)s
  shiftsol = Psi0;
  shiftsol(2:end,:) = shiftsol(1:end-1,:);

  H2mu = ttt(tensor(ones(nx,1),nx),tensor(mu/dx,nmu));
  H2 = H2mu - sigma/2;

  rhs = H2.*shiftsol + Q;
  rhs(1,1) = 0;
  rhs(nx,nmu) = 0;

  %convert to linear system and solve
  Hmat = reshape(permute(H1,[1,3,2,4]),[nx*nmu,nx*nmu]);
  Psivec = double(Hmat)\double(rhs(:));
  Psi1 = tensor(reshape(Psivec,[nx,nmu]));
  converr = norm(Psi1-Psi0)/norm(Psi0);
  fprintf('iter = %3.d, converr = %.5e \n', iter, converr);
  if converr<=convtol
    break;
  else
    Psi0 = Psi1;
%     disp(Psi0)
  end
end




