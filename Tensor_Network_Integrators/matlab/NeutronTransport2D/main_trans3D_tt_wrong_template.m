close all; clear; clc;
run setup2DSlab.m

%% Slab Dimensions
x0 = 0.0; x1 = 4.2;
y0 = 0.0; y1 = 4.2;
z0 = 0.0; z1 = 4.2;
Param = load('./Data/OneGrpData.mat');
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
nz = 3; %z
L = 4; %total number of points is 2^d*(L/2)^2
% nmu = 2; % mu
% neta = 2; % eta
% nxi = 2; % xi

dx = (x1-x0)/(nx-1); %spacial step
dy = (y1-y0)/(ny-1); %spacial step
dz = (z1-z0)/(nz-1); %spacial step

% % regenerate mu and weight
idimen = 3;
octs   = 2^idimen;
nords = octs*(L/2)^2;
% [ mu, eta, xi, wgt ] = GetSNOrdinates( L, idimen);
[ mu, eta, xi, wgt ] = GetChebyLegendre ( L, idimen );
%% forming components
%positive mu and eta
difmatp = 1/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));
% angp = diag([zeros(1,nmu/2), mu(nmu/2+1:nmu)]);
intpmatp = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1));
intgmat = ones(numel(mu),1)*wgt;

intpmatp_noBC = intpmatp;
intpmatp_noBC(1,:) = 0;

%negative mu and eta
difmatm = 1/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
% angm = diag([mu(1:nmu/2), zeros(1,nmu/2)]);
intpmatm = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),1));
% intgm = flipud(intgp);

intpmatm_noBC = intpmatm;
intpmatm_noBC(end,:) = 0;

% collect
difmat ={difmatm,difmatp};
% ang = {angm, angp};
intp = {intpmatm,intpmatp};
% intg = {intgm,intgp};

intp_noBC = {intpmatm_noBC, intpmatp_noBC};

%% Now build the tensor
ttH = 0; %LHS operator
ttS = 0; %RHS scattering operator
ttF = 0; % RHS fission operator

for imu = 1:2 % 1-minus, 2-plus
  for jeta = 1:2
    for kxi = 1:2
      %d/dx term
      Ldx = matrices_to_tt_matrix_fn({difmat{imu}, intp{jeta}, intp{kxi}, ...
        ang{imu},ang{jeta}~=0, ang{kxi}~=0, eye(nE)});

      %d/dy term
      Ldy = matrices_to_tt_matrix_fn({intp{imu}, difmat{jeta},intp{kxi},...
        ang{imu}~=0, ang{jeta}, ang{kxi}~=0, eye(nE)});

      %d/dz term
      Ldz = matrices_to_tt_matrix_fn({intp{imu}, intp{jeta}, difmat{kxi},...
        ang{imu}~=0, ang{jeta}~=0, ang{kxi}, eye(nE)});
      
      %interpolation term
      Lintp = matrices_to_tt_matrix_fn({intp{imu}, intp{jeta}, intp{kxi},...
        ang{imu}~=0, ang{jeta}~=0, ang{kxi}~=0, diag(sigt)});
      
      Htemp = Ldx + Ldy + Ldz + Lintp;
      
      %scattering
      Stemp = matrices_to_tt_matrix_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        intg{imu}, intg{jeta}, intg{kxi}, sigs});
      %fission
      Ftemp = matrices_to_tt_matrix_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        intg{imu}, intg{jeta}, intg{kxi}, nusigf});

      %accumulate
      if ~isa(ttH,'tt_matrix')
        ttH = Htemp;
        ttS = Stemp;
        ttF = Ftemp;
      else
        ttH = ttH + Htemp;
        ttS = ttS + Stemp;
        ttF = ttF + Ftemp;
      end
    end
  end
end
%% Rounding operators
tt_tol = 1e-6;
ttH = round(ttH, tt_tol);
ttHmat = full(ttH);
ttS = round(ttS, tt_tol);
ttF = round(ttF, tt_tol);
%% solve eigenvalue in tt
fixed_point_tol = 1e-6;

niter = 100;
[ktt, ttPsi1] = tt_fixed_point_eig_solve(ttH, ttS, ttF, niter,...
  fixed_point_tol, tt_tol);

fprintf('eigenvalue k = %.10f \n', ktt)
fprintf('compression ratio of H = %.5e \n', compress_ratio_tt(ttH));
fprintf('compression ratio of S = %.5e \n', compress_ratio_tt(ttS));
fprintf('compression ratio of F = %.5e \n', compress_ratio_tt(ttF));









