close all; clear; clc;
if ~isdeployed
  run setup2DSlab.m
end

%% Slab Dimensions
x0 = 0.0; x1 = 10;
y0 = 0.0; y1 = 10;
z0 = 0.0; z1 = 10;
Param = load('NTE_3D.mat');

%%
chi = Param.chi;
nusigf = Param.nusigf;
sigt = Param.sigt;
sigs = Param.sigs;

nusigf = chi.*nusigf';
nE = length( sigt );
fixed_point_tol = 1e-3;
tt_tol = fixed_point_tol*0.01;

%% ASSUME Nx=Ny, dx=dy
n = 128;
nx = n; %x
ny = n; %y
nz = n; %z
L = 8;

dx = (x1-x0)/(nx-1); %spacial step
dy = (y1-y0)/(ny-1); %spacial step
dz = (z1-z0)/(nz-1); %spacial step

% % regenerate mu and weight
idimen = 3;
octs   = 2^idimen;
nords = octs*(L/2)^2;
[ mu, eta, xi, wgt ] = GetChebyLegendre ( L, idimen );

fprintf('Problem grid size = %d \n',n)

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
      Ldx = matrices_to_qtt_matrix_fn({difmat{imu}, intp{jeta}, intp{kxi}, ...
        angmat, eye(nE)},tt_tol);
      
      % %d/dy term
      angmat = kron(Itemp,(jeta-1.5)*2*diag(eta));
      Ldy = matrices_to_qtt_matrix_fn({intp{imu}, difmat{jeta},intp{kxi},...
        angmat, eye(nE)},tt_tol);
      %d/dz term
      angmat = kron(Itemp,(kxi-1.5)*2*diag(xi));
      Ldz = matrices_to_qtt_matrix_fn({intp{imu}, intp{jeta}, difmat{kxi},...
        angmat, eye(nE)},tt_tol);

      %interpolation term
      angmat = kron(Itemp, eye(numel(mu)));
      Lintp = matrices_to_qtt_matrix_fn({intp{imu}, intp{jeta}, intp{kxi},...
        angmat, diag(sigt)},tt_tol);

      Htemp = Ldx + Ldy + Ldz + Lintp;
      
      %scattering
      %construct the operator for the integral
      Itemp(octidx,:) = 1;
      Stemp = matrices_to_qtt_matrix_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        kron(Itemp,intgmat), sigs},tt_tol);

      %fission
      Ftemp = matrices_to_qtt_matrix_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        kron(Itemp,intgmat), nusigf},tt_tol);

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
ttH = round(ttH, tt_tol);
ttS = round(ttS, tt_tol);
ttF = round(ttF, tt_tol);
%% solve eigenvalue in tt
niter = 100;
mv_nswp = 20;
amen_solve_nswp = 10;
kickrank = 2;
params = struct('niter',niter,'epsi',fixed_point_tol,'tt_tol',tt_tol);
mv_opts = {'verb',0,'nswp',mv_nswp};
amen_solve_opts = {'verb',0,'nswp', amen_solve_nswp,'kickrank',kickrank,...
  'trunc_norm','fro','rmax', 120, 'x0',0};

[ktt, ttPsi1, Itertime, convrate,lambda] = ...
  tt_fixed_point_eig_solve(ttH, ttS, ttF, params, mv_opts, amen_solve_opts); 
%% PRINT THE RESULTS
true_k = 1.02758273;
fprintf('eigenvalue k = %.10f \n', ktt)
fprintf('eigenvalue errors = %.5e \n', norm(ktt-true_k));
fprintf('compression ratio of H = %.5e \n', compress_ratio_tt(ttH));
fprintf('compression ratio of S = %.5e \n', compress_ratio_tt(ttS));
fprintf('compression ratio of F = %.5e \n', compress_ratio_tt(ttF));
fprintf('Elapsed time = %.2f seconds \n',sum(Itertime))

%%
filename = sprintf('./Results/keffective_darwin_%d.mat',n);
save(filename, ...
  'ktt','ttPsi1','ttH','ttS','ttF','Itertime','true_k','params','mv_opts',...
  'amen_solve_opts','convrate','lambda');

