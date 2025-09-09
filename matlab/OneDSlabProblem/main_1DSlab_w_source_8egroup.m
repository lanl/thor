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


% RHS1 Q
Q = 1*ones(nx,nmu,nE); %this is a vector of ones
Q(1,1,:) = 0; %BC mu pos
Q(nx,nmu,:) = 0; % BC mu neg

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
%% convert to linear system and solve
Hmat = reshape(permute(H,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]);
Qvec = Q(:);

Psivec = double(Hmat)\double(Qvec);
Psi = reshape(Psivec,[nx,nmu,nE]);

%% solve by tensor train
H_for_tt = permute(double(H),[1,3,5,2,4,6]);
Htt = tt_matrix(double(H_for_tt));

Httfull = reshape(full(Htt),size(H_for_tt));
temperr = norm(tensor(H_for_tt-Httfull))/norm(tensor(H_for_tt));

fprintf('Error between Htt and H = %.5e \n', temperr)
fprintf('compress ratio in H = %.2f \n',compress_ratio_tt(Htt,Htt.r, prod(size(Htt),2)));
Qtt = tt_tensor(Q);

Psitt = amen_solve2(Htt, Qtt, 1e-6);
solerr = norm(Psi(:) - full(Psitt));

fprintf('Difference between Psi and Psitt = %.5e \n', solerr);
