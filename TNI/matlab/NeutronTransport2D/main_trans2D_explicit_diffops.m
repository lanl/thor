close all; clear; clc;
run setup2DSlab.m

%% Slab Dimensions
x0 = 0.0; x1 = 4.2;
y0 = 0.0; y1 = 4.2;

% Total Cross Section
sigt = 0.45391;
% sigt = 0;
% sigt   = [ 2.3336e-01 2.3105e-01 3.0084e-01 2.8809e-01 3.8538e-01 5.0120e-01 6.1661e-01 7.8866e-01];
%sigt   = [ 2.9138e-1 4.7284e-1 ];
% sigs = 0.34612;
% nusigf = 0.24683;

% mu   = [ -0.57735027, 0.57735027, -0.57735027, 0.57735027 ];
mu = [-0.57735027, 0.57735027];
% eta  = [ -0.57735027, -0.57735027, 0.57735027, 0.57735027 ];
eta = [-0.57735027, 0.57735027];

% wgt  = [ 0.25, 0.25, 0.25, 0.25 ];
wgt = [0.5, 0.5];

% Problem Parameters
%% ASSUME Nx=Ny, dx=dy
nx = 3; %x
ny = 3; %y
nmu = 2; % mu
neta = 2; % eta
nE = length( sigt );

dx = (x1-x0)/(nx-1); %spacial step
dy = (x1-x0)/(ny-1); %spacial step


%% forming H
%positive mu and eta
difmatp = 1/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));
angp = diag([zeros(1,nmu/2), mu(nmu/2+1:nmu)]);
intpmatp = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1));

%negative mu and eta
difmatm = 1/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
angm = diag([mu(1:nmu/2), zeros(1,nmu/2)]);
intpmatm = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),1));


%% testing the differential part
dxten = ttt(tensor(difmatm),tensor(intpmatm));
dxmat = reshape(permute(dxten,[1,3,2,4]),[nx*ny,nx*ny]);

dyten = ttt(tensor(intpmatm), tensor(difmatm));
dymat = reshape(permute(dyten,[1,3,2,4]),[nx*ny,nx*ny]);

%% (1) mu<0, eta<0
Ldx = ttt_out_prod_fn({difmatm,intpmatm,angm,angm~=0});
Ldy = ttt_out_prod_fn({intpmatm, difmatm, angm~=0, angm});

Lmm = Ldx + Ldy;
Lmmmat = get_mat_8D_fn(Lmm,nx,ny,nmu,neta);

%% (2) mu>0 (x>0), eta<0 (y<0)
Ldx = ttt_out_prod_fn({difmatp,intpmatm,angp,angm~=0});
Ldy = ttt_out_prod_fn({intpmatp, difmatm, angp~=0, angm});

Lpm = Ldx + Ldy;
Lpmmat = get_mat_8D_fn(Lpm,nx,ny,nmu,neta);

%% (3) mu<0 (x<0), eta>0 (y>0)
Ldx = ttt_out_prod_fn({difmatm,intpmatp,angm,angp~=0});
Ldy = ttt_out_prod_fn({intpmatm, difmatp, angm~=0, angp});

Lmp = Ldx + Ldy;
Lmpmat = get_mat_8D_fn(Lmp,nx,ny,nmu,neta);

%% (4) mu>0 (x>0), eta>0 (y>0)
Ldx = ttt_out_prod_fn({difmatp,intpmatp,angp,angp~=0});
Ldy = ttt_out_prod_fn({intpmatp, difmatp, angp~=0, angp});

Lpp = Ldx + Ldy;
Lppmat = get_mat_8D_fn(Lpp,nx,ny,nmu,neta);
%%
H = Lpp + Lpm + Lmp + Lmm;
Hmat = get_mat_8D_fn(H,nx,ny,nmu,neta);




