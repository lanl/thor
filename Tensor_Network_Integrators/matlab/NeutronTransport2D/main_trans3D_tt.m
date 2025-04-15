close all; clear; clc;
run setup2DSlab.m

%{ 
In this code, the quadrature rules for ordinates are on a sphere
All 3 angles will be merged into one long dimension
Psi is a 5D tensor and have the dimensions: nx x ny x nz x nords x nG
%}
%% Slab Dimensions
x0 = 0.0; x1 = 4.2;
y0 = 0.0; y1 = 4.2;
z0 = 0.0; z1 = 4.2;

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

%% 

nx = 3; %x
ny = 3; %y
nz = 3; %z
L = 4; %total number of points is 2^d*(L/2)^2

dx = (x1-x0)/(nx-1); %spacial step
dy = (y1-y0)/(ny-1); %spacial step
dz = (z1-z0)/(nz-1); %spacial step

% % regenerate mu and weight
idimen = 3;
octs   = 2^idimen;
nords = octs*(L/2)^2;
[ mu, eta, xi, wgt ] = GetChebyLegendre ( L, idimen );

%% Source Q -- deal later
% Q = 1*ones(nx,ny,nmu,neta,nE);
% % 4 boundary conditions
% Q(nx,:,1,:,:) = 0; %mm
% Q(:,ny,:,1,:) = 0; %mm
% Q(1,:,nmu,:,:) = 0; %pm
% Q(:,ny,:,1,:) = 0; %pm
% Q(:,1,:,neta,:) = 0; %mp
% Q(nx,:,1,:,:) = 0; %mp
% Q(:,1,:,neta,:) = 0; %pp
% Q(1,:,nmu,:,:) = 0; %pp

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

for kxi = 1:2 % 1-minus, 2-plus
  for jeta = 1:2
    for imu = 1:2
      Itemp = zeros(octs,octs);
      octidx = 4*(kxi-1) + 2*(jeta-1) + imu;
      % fprintf('octant %d \n',octidx)
      Itemp(octidx,octidx) = 1 ; %map imu from [1,2] to [-1,1]

      %d/dx term
      angmat = kron(Itemp,(imu-1.5)*2*diag(mu));
      Ldx = matrices_to_tt_matrix_fn({difmat{imu}, intp{jeta}, intp{kxi}, ...
        angmat, eye(nE)});
      
      % %d/dy term
      angmat = kron(Itemp,(jeta-1.5)*2*diag(eta));
      Ldy = matrices_to_tt_matrix_fn({intp{imu}, difmat{jeta},intp{kxi},...
        angmat, eye(nE)});
      %d/dz term
      angmat = kron(Itemp,(kxi-1.5)*2*diag(xi));
      Ldz = matrices_to_tt_matrix_fn({intp{imu}, intp{jeta}, difmat{kxi},...
        angmat, eye(nE)});

      %interpolation term
      angmat = kron(Itemp, eye(numel(mu)));
      Lintp = matrices_to_tt_matrix_fn({intp{imu}, intp{jeta}, intp{kxi},...
        angmat, diag(sigt)});

      Htemp = Ldx + Ldy + Ldz + Lintp;
      
      %scattering
      %construct the operator for the integral
      Itemp(octidx,:) = 1;
      Stemp = matrices_to_tt_matrix_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        kron(Itemp,intgmat), sigs});

      %fission
      Ftemp = matrices_to_tt_matrix_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        kron(Itemp,intgmat), nusigf});

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

%%
%check error with ref data
if 0
  refdata = load('./References/H_cartesian_2angles.mat');
  Hmatref = refdata.Hmat;
  Smatref = refdata.Smat;
  Fmatref = refdata.Fmat;

  norm(Hmat-Hmatref)
  norm(Smat-Smatref)
  norm(Hmat-Hmatref)
end

%% Rounding operators
tt_tol = 1e-6;
ttH = round(ttH, tt_tol);
ttS = round(ttS, tt_tol);
ttF = round(ttF, tt_tol);
%% solve eigenvalue in tt
fixed_point_tol = 1e-6;

niter = 100;
[ktt, ttPsi1, Itertime, convrate, lambda] = tt_fixed_point_eig_solve(ttH, ttS, ttF, niter,...
  fixed_point_tol, tt_tol);

fprintf('eigenvalue k = %.10f \n', ktt)
fprintf('compression ratio of H = %.5e \n', compress_ratio_tt(ttH));
fprintf('compression ratio of S = %.5e \n', compress_ratio_tt(ttS));
fprintf('compression ratio of F = %.5e \n', compress_ratio_tt(ttF));








