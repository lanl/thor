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

%Param = load('./Data/OneGrpData.mat');
Param = load('./Data/EightGrpData.mat');

sigt = Param.sigt;
sigs = Param.sigs;
nusigf = 3*Param.nusigf;
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
L = 2; %total number of points is 2^d*(L/2)^2
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
H = 0; %LHS operator
S = 0; %RHS scattering operator
F = 0; % RHS fission operator

for kxi = 1:2 % 1-minus, 2-plus
  for jeta = 1:2
    for imu = 1:2

      Itemp = zeros(octs,octs);
      octidx = 4*(kxi-1) + 2*(jeta-1) + imu;
      % fprintf('octant %d \n',octidx)
      Itemp(octidx,octidx) = 1 ; %map imu from [1,2] to [-1,1]

      %d/dx term
      angmat = kron(Itemp,(imu-1.5)*2*diag(mu));
      Ldx = ttt_out_prod_fn({difmat{imu}, intp{jeta}, intp{kxi}, ...
        angmat, eye(nE)});
      Ldxmat = get_mat_10D_fn(Ldx,nx,ny,nz,nords,nE);
      
      % %d/dy term
      angmat = kron(Itemp,(jeta-1.5)*2*diag(eta));
      Ldy = ttt_out_prod_fn({intp{imu}, difmat{jeta},intp{kxi},...
        angmat, eye(nE)});
      Ldymat = get_mat_10D_fn(Ldy,nx,ny,nz,nords,nE);
      %d/dz term
      angmat = kron(Itemp,(kxi-1.5)*2*diag(xi));
      Ldz = ttt_out_prod_fn({intp{imu}, intp{jeta}, difmat{kxi},...
        angmat, eye(nE)});
      Ldzmat = get_mat_10D_fn(Ldz,nx,ny,nz,nords,nE);

      %interpolation term
      angmat = kron(Itemp, eye(numel(mu)));
      Lintp = ttt_out_prod_fn({intp{imu}, intp{jeta}, intp{kxi},...
        angmat, diag(sigt)});

      Htemp = Ldx + Ldy + Ldz + Lintp;
      Htempmat = get_mat_10D_fn(Htemp,nx,ny,nz,nords,nE);
      

      %scattering
      %construct the operator for the integral
      Itemp(octidx,:) = 1;
      Stemp = ttt_out_prod_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        kron(Itemp,intgmat), sigs});

      %fission
      Ftemp = ttt_out_prod_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        kron(Itemp,intgmat), nusigf});

      %accumulate
      H = H + Htemp;
      S = S + Stemp;
      F = F + Ftemp;
    end
  end
end

%% solve linear system
Hmat = get_mat_10D_fn(H,nx,ny,nz,nords,nE);
Smat = get_mat_10D_fn(S,nx,ny,nz,nords,nE);
Fmat = get_mat_10D_fn(F,nx,ny,nz,nords,nE);

%% 
% %check error with ref data
% if 0
%   refdata = load('./References/H_cartesian_2angles.mat');
%   Hmatref = refdata.Hmat;
%   Smatref = refdata.Smat;
%   Fmatref = refdata.Fmat;
% 
%   norm(Hmat-Hmatref)
%   norm(Smat-Smatref)
%   norm(Hmat-Hmatref)
% end
%%
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
%%
%% Rounding operators
tt_tol = 1e-6;
ttH = round(ttH, tt_tol);
ttHmat = full(ttH);
ttS = round(ttS, tt_tol);
ttF = round(ttF, tt_tol);
%%
niter = 100;
mv_nswp = 20;
amen_solve_nswp = 10;
kickrank = 2;
params = struct('niter',niter,'epsi',fixed_point_tol,'tt_tol',tt_tol);
mv_opts = {'verb',0,'nswp',mv_nswp};
amen_solve_opts = {'verb',0,'nswp', amen_solve_nswp,'kickrank',kickrank,...
  'trunc_norm','fro','rmax', 120, 'x0',0};
[ktt, ttPsi1] = tt_fixed_point_eig_solve(ttH, ttS, ttF, ...
  params, mv_opts, amen_solve_opts);

%%
fprintf('TT eigenvalue = %.10f \n',ktt)
fprintf('Error in eigenvalue = %.2e \n',abs(ktt-k))
fprintf('compression ratio of H = %.5e \n', compress_ratio_tt(ttH));
fprintf('compression ratio of S = %.5e \n', compress_ratio_tt(ttS));
fprintf('compression ratio of F = %.5e \n', compress_ratio_tt(ttF));








