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
b = 4.2;
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
c = mu/dx  + sigma/2;
d = -mu/dx + sigma/2;

Lpos = diag(c(1)*ones(nx,1)) + diag(d(1:nmu/2)*ones(nx-1,1),-1);
Lpos(1,1) = 1;
Lneg = diag(d(2)*ones(nx,1)) + diag(c(nmu/2+1:end)*ones(nx-1,1),1);
Lneg(nx,nx) = 1;
Lpos = tensor(Lpos);
Lneg = tensor(Lneg);

H1 = ttt(tensor(Lpos), tensor([1,0;0,0]));
H2 = ttt(tensor(Lneg), tensor([0,0;0,1]));
H = H1+H2;

% adjust the boundary with coefficients c and d
Q(1,1) = Q(1,1)/c(1);
Q(nx,nmu) = Q(nx,nmu)/d(nmu);

%% check the construction of H
mu_mat = diag(mu)/dx;

Hmu = ttt(tensor(eye(nx)),tensor(mu_mat));
Hg = sigma/2.*ttt(tensor(eye(nx)), tensor(eye(nmu)));
L = ttt(tensor(diag(ones(nx-1,1),-1)),tensor(eye(nmu)));
Hmu1 = ttt(tensor(diag(ones(nx-1,1),+1)),tensor(mu_mat));
Hg1 = sigma/2.*ttt(tensor(diag(ones(nx-1,1),+1)), tensor(eye(nmu)));
Htemp = Hmu + Hg - Hg1 + Hmu1;
%% convert to linear system and solve
Hmat = reshape(permute(H,[1,3,2,4]),[nx*nmu,nx*nmu]);
Qvec = Q(:);

Psivec = double(Hmat)\double(Qvec);
Psi = reshape(Psivec,[nx,nmu])
