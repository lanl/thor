close all; clear; clc;
run setup2DSlab.m

%% Slab Dimensions
x0 = 0.0; x1 = 4.2;
y0 = 0.0; y1 = 4.2;
% Param = load('./Data/OneGrpData.mat');
Param = load('./Data/EightGrpData.mat');

sigt = Param.sigt;
sigs = Param.sigs;
nusigf = Param.nusigf;
nE = length( sigt );
chi    = ones( nE,1 )/nE;
nusigf = chi.*nusigf';

% mu = [-0.57735027, 0.57735027];
% eta = [-0.57735027, 0.57735027];
% wgt = [0.5, 0.5];

%% ASSUME Nx=Ny, dx=dy
nx = 3; %x
ny = 3; %y
nmu = 2; % mu
neta = 2; % eta


dx = (x1-x0)/(nx-1); %spacial step
dy = (x1-x0)/(ny-1); %spacial step

% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]); %assuming eta is the same as mu here
wgt = flip(wgt'/2);

%% Source Q%
Q = 1*ones(nx,ny,nmu,neta,nE);
% 4 boundary conditions
Q(nx,:,1,:,:) = 0; %mm
Q(:,ny,:,1,:) = 0; %mm
Q(1,:,nmu,:,:) = 0; %pm
Q(:,ny,:,1,:) = 0; %pm
Q(:,1,:,neta,:) = 0; %mp
Q(nx,:,1,:,:) = 0; %mp
Q(:,1,:,neta,:) = 0; %pp
Q(1,:,nmu,:,:) = 0; %pp

%% forming components
%positive mu and eta
difmatp = 1/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));
angp = diag([zeros(1,nmu/2), mu(nmu/2+1:nmu)]);
intpmatp = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1));
intgp = [zeros(nmu/2,1);ones(nmu/2,1)]*wgt';

intpmatp_noBC = intpmatp;
intpmatp_noBC(1,:) = 0;

%negative mu and eta
difmatm = 1/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
angm = diag([mu(1:nmu/2), zeros(1,nmu/2)]);
intpmatm = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),1));
intgm = flipud(intgp);

intpmatm_noBC = intpmatm;
intpmatm_noBC(end,:) = 0;

%% (1) mu<0, eta<0
%d/dx term
Ldx = ttt_out_prod_fn({difmatm,intpmatm,angm,angm~=0,eye(nE)});
%d/dy term
Ldy = ttt_out_prod_fn({intpmatm, difmatm, angm~=0, angm, eye(nE)});
%interpolation term
Lintp = ttt_out_prod_fn({intpmatm,intpmatm,angm~=0,angm~=0, diag(sigt)});

Lmm = Ldx + Ldy + Lintp;
Lmmmat = get_mat_10D_fn(Lmm,nx,ny,nmu,neta,nE);
disp(Lmmmat(1:9,1:9));

%scattering
Smm = ttt_out_prod_fn({intpmatm_noBC,intpmatm_noBC,intgm,intgm,sigs});
% Smmmat = get_mat_10D_fn(Smm,nx,ny,nmu,neta,nE);

%fission
Fmm = ttt_out_prod_fn({intpmatm_noBC,intpmatm_noBC,intgm,intgm,nusigf});
% Fmmmat = get_mat_10D_fn(Fmm,nx,ny,nmu,neta,nE);

%% (2) mu>0 (x>0), eta<0 (y<0)
%d/dx term
Ldx = ttt_out_prod_fn({difmatp,intpmatm,angp,angm~=0,eye(nE)});
%d/dy term
Ldy = ttt_out_prod_fn({intpmatp, difmatm, angp~=0, angm,eye(nE)});
%interpolation term
Lintp = ttt_out_prod_fn({intpmatp,intpmatm,angp~=0,angm~=0, diag(sigt)});

Lpm = Ldx + Ldy + Lintp; %for H
Lpmmat = get_mat_10D_fn(Lpm,nx,ny,nmu,neta,nE);

%scattering
Spm = ttt_out_prod_fn({intpmatp_noBC,intpmatm_noBC,intgp,intgm,sigs});

%fission
Fpm = ttt_out_prod_fn({intpmatp_noBC,intpmatm_noBC,intgp,intgm,nusigf});



%% (3) mu<0 (x<0), eta>0 (y>0)
%d/dx term
Ldx = ttt_out_prod_fn({difmatm,intpmatp,angm,angp~=0,eye(nE)});
%d/dy term
Ldy = ttt_out_prod_fn({intpmatm, difmatp, angm~=0, angp,eye(nE)});
%interpolation term
Lintp = ttt_out_prod_fn({intpmatm,intpmatp,angm~=0,angp~=0, diag(sigt)});

Lmp = Ldx + Ldy + Lintp;
Lmpmat = get_mat_10D_fn(Lmp,nx,ny,nmu,neta,nE);

%scattering
Smp = ttt_out_prod_fn({intpmatm_noBC,intpmatp_noBC,intgm,intgp,sigs});

%fission
Fmp = ttt_out_prod_fn({intpmatm_noBC,intpmatp_noBC,intgm,intgp,nusigf});

%% (4) mu>0 (x>0), eta>0 (y>0)
%d/dx term
Ldx = ttt_out_prod_fn({difmatp,intpmatp,angp,angp~=0,eye(nE)});
%d/dy term
Ldy = ttt_out_prod_fn({intpmatp, difmatp, angp~=0, angp,eye(nE)});
%interpolation term
Lintp = ttt_out_prod_fn({intpmatp,intpmatp,angp~=0,angp~=0, diag(sigt)});

Lpp = Ldx + Ldy+Lintp;
Lppmat = get_mat_10D_fn(Lpp,nx,ny,nmu,neta,nE);

%scattering
Spp = ttt_out_prod_fn({intpmatp_noBC,intpmatp_noBC,intgp,intgp,sigs});

%fission
Fpp = ttt_out_prod_fn({intpmatp_noBC,intpmatp_noBC,intgp,intgp,nusigf});
%% Building H
H = Lpp + Lpm + Lmp + Lmm;
S = Spp + Spm + Smp + Smm;
F = Fpp + Fpm + Fmp + Fmm;

% H_for_tt = permute(double(H),[1,3,5,7,9,2,4,6,8,10]);
% Htt = tt_matrix(H_for_tt);
% Qtt = tt_tensor(Q);
% Psitt = amen_solve2(Htt,Qtt,1e-6)
%% solve linear system
Hmat = get_mat_10D_fn(H,nx,ny,nmu,neta,nE);
Smat = get_mat_10D_fn(S,nx,ny,nmu,neta,nE);
Fmat = get_mat_10D_fn(F,nx,ny,nmu,neta,nE);
% Hmat = get_mat_10D_fn(H,nx,ny,nmu,neta,nE);
% Qvec = Q(:);
% Psivec = (Hmat-Smat-Fmat)\Qvec;
%% solve eigenvalue problem
A = (Hmat-Smat)\Fmat;
k = max(abs(eig(A)));

%% solve eigenvalue in tt
ttH = convert_tensor_op_to_tt_matrix_fn(H,10);
ttS = convert_tensor_op_to_tt_matrix_fn(S,10);
ttF = convert_tensor_op_to_tt_matrix_fn(F,10);
fixed_point_tol = 1e-6;
tol = 1e-6;
niter = 100;
[ktt, ttPsi1] = tt_fixed_point_eig_solve(ttH, ttS, ttF, niter,...
  fixed_point_tol,tol);

fprintf('TT eigenvalue = %.10f \n',ktt)
fprintf('Error in eigenvalue = %.2e \n',abs(ktt-k))








