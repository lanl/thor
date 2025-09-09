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

tt_tol = 1e-6;

%% ASSUME Nx=Ny, dx=dy
nx = 8; %x
ny = 8; %y
nmu = 8; % mu
neta = 8; % eta

dx = (x1-x0)/(nx-1); %spacial step
dy = (x1-x0)/(ny-1); %spacial step

% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]); %assuming eta is the same as mu here
wgt = flip(wgt'/2);
mu = mu';

%% Source Q
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

Ldx = matrices_to_qtt_matrix_fn({difmatm,intpmatm,angm,angm~=0,eye(nE)}, tt_tol);
%d/dy term
Ldy = matrices_to_qtt_matrix_fn({intpmatm, difmatm, angm~=0, angm, eye(nE)}, tt_tol);
%interpolation term
Lintp = matrices_to_qtt_matrix_fn({intpmatm,intpmatm,angm~=0,angm~=0, diag(sigt)},tt_tol);

Lmm = Ldx + Ldy + Lintp;

%scattering
Smm = matrices_to_qtt_matrix_fn({intpmatm_noBC,intpmatm_noBC,intgm,intgm,sigs},tt_tol);
% Smmmat = get_mat_10D_fn(Smm,nx,ny,nmu,neta,nE);

%fission
Fmm = matrices_to_qtt_matrix_fn({intpmatm_noBC,intpmatm_noBC,intgm,intgm,nusigf},tt_tol);
% Fmmmat = get_mat_10D_fn(Fmm,nx,ny,nmu,neta,nE);

%% (2) mu>0 (x>0), eta<0 (y<0)
%d/dx term
Ldx = matrices_to_qtt_matrix_fn({difmatp,intpmatm,angp,angm~=0,eye(nE)}, tt_tol);
%d/dy term
Ldy = matrices_to_qtt_matrix_fn({intpmatp, difmatm, angp~=0, angm,eye(nE)}, tt_tol);
%interpolation term
Lintp = matrices_to_qtt_matrix_fn({intpmatp,intpmatm,angp~=0,angm~=0, diag(sigt)},tt_tol);

Lpm = Ldx + Ldy + Lintp; %for H

%scattering
Spm = matrices_to_qtt_matrix_fn({intpmatp_noBC,intpmatm_noBC,intgp,intgm,sigs},tt_tol);

%fission
Fpm = matrices_to_qtt_matrix_fn({intpmatp_noBC,intpmatm_noBC,intgp,intgm,nusigf},tt_tol);



%% (3) mu<0 (x<0), eta>0 (y>0)
%d/dx term
Ldx = matrices_to_qtt_matrix_fn({difmatm,intpmatp,angm,angp~=0,eye(nE)}, tt_tol);
%d/dy term
Ldy = matrices_to_qtt_matrix_fn({intpmatm, difmatp, angm~=0, angp,eye(nE)}, tt_tol);
%interpolation term
Lintp = matrices_to_qtt_matrix_fn({intpmatm,intpmatp,angm~=0,angp~=0, diag(sigt)}, tt_tol);

Lmp = Ldx + Ldy + Lintp;

%scattering
Smp = matrices_to_qtt_matrix_fn({intpmatm_noBC,intpmatp_noBC,intgm,intgp,sigs}, tt_tol);

%fission
Fmp = matrices_to_qtt_matrix_fn({intpmatm_noBC,intpmatp_noBC,intgm,intgp,nusigf}, tt_tol);

%% (4) mu>0 (x>0), eta>0 (y>0)
%d/dx term
Ldx = matrices_to_qtt_matrix_fn({difmatp,intpmatp,angp,angp~=0,eye(nE)}, tt_tol);
%d/dy term
Ldy = matrices_to_qtt_matrix_fn({intpmatp, difmatp, angp~=0, angp,eye(nE)}, tt_tol);
%interpolation term
Lintp = matrices_to_qtt_matrix_fn({intpmatp,intpmatp,angp~=0,angp~=0, diag(sigt)}, tt_tol);

Lpp = Ldx + Ldy+Lintp;

%scattering
Spp = matrices_to_qtt_matrix_fn({intpmatp_noBC,intpmatp_noBC,intgp,intgp,sigs}, tt_tol);

%fission
Fpp = matrices_to_qtt_matrix_fn({intpmatp_noBC,intpmatp_noBC,intgp,intgp,nusigf}, tt_tol);
%% Building H
ttH = Lpp + Lpm + Lmp + Lmm;
ttS = Spp + Spm + Smp + Smm;
ttF = Fpp + Fpm + Fmp + Fmm;

ttH = round(ttH, tt_tol);
ttS = round(ttS, tt_tol);
ttF = round(ttF, tt_tol);

fprintf('compression ratio of H = %.5e \n', compress_ratio_tt(ttH));
fprintf('compression ratio of S = %.5e \n', compress_ratio_tt(ttS));
fprintf('compression ratio of F = %.5e \n', compress_ratio_tt(ttF));
%% solve linear system
% Hmat = get_mat_10D_fn(ttH,nx,ny,nmu,neta,nE);
% Smat = get_mat_10D_fn(ttS,nx,ny,nmu,neta,nE);
% Fmat = get_mat_10D_fn(ttF,nx,ny,nmu,neta,nE);
% 
% %% solve eigenvalue problem
% A = (Hmat-Smat)\Fmat;
% k = max(abs(eig(A)));

%% solve eigenvalue in tt
fixed_point_tol = 1e-6;

niter = 100;
[ktt, ttPsi1] = tt_fixed_point_eig_solve(ttH, ttS, ttF, niter,...
  fixed_point_tol, tt_tol);

fprintf('eigenvalue k = %.10f \n', ktt)










