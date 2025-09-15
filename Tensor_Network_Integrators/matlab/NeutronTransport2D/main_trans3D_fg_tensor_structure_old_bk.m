close all; clear; clc;
run setup2DSlab.m

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
nmu = 2; % mu
neta = 2; % eta
nxi = 2; % xi

dx = (x1-x0)/(nx-1); %spacial step
dy = (y1-y0)/(ny-1); %spacial step
dz = (z1-z0)/(nz-1); %spacial step

% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]); %assuming eta is the same as mu here
wgt = flip(wgt'/2);

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

% collect
difmat ={difmatm,difmatp};
ang = {angm, angp};
intp = {intpmatm,intpmatp};
intg = {intgm,intgp};

intp_noBC = {intpmatm_noBC, intpmatp_noBC};

%% Now build the tensor
H = 0; %LHS operator
S = 0; %RHS scattering operator
F = 0; % RHS fission operator

for imu = 1:2 % 1-minus, 2-plus
  for jeta = 1:2
    for kxi = 1:2
      %d/dx term
      Ldx = ttt_out_prod_fn({difmat{imu}, intp{jeta}, intp{kxi}, ...
        ang{imu},ang{jeta}~=0, ang{kxi}~=0, eye(nE)});

      %d/dy term
      Ldy = ttt_out_prod_fn({intp{imu}, difmat{jeta},intp{kxi},...
        ang{imu}~=0, ang{jeta}, ang{kxi}~=0, eye(nE)});

      %d/dz term
      Ldz = ttt_out_prod_fn({intp{imu}, intp{jeta}, difmat{kxi},...
        ang{imu}~=0, ang{jeta}~=0, ang{kxi}, eye(nE)});
      
      %interpolation term
      Lintp = ttt_out_prod_fn({intp{imu}, intp{jeta}, intp{kxi},...
        ang{imu}~=0, ang{jeta}~=0, ang{kxi}~=0, diag(sigt)});
      
      Htemp = Ldx + Ldy + Ldz + Lintp;
      
      %scattering
      Stemp = ttt_out_prod_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        intg{imu}, intg{jeta}, intg{kxi}, sigs});
      %fission
      Ftemp = ttt_out_prod_fn({intp_noBC{imu}, intp_noBC{jeta}, intp_noBC{kxi},...
        intg{imu}, intg{jeta}, intg{kxi}, nusigf});
      kron(kron(intg{imu},intg{jeta}),intg{kxi})
      %accumulate
      H = H + Htemp;
      S = S + Stemp;
      F = F + Ftemp;
    end
  end
end

%% solve linear system
Hmat = get_mat_14D_fn(H,nx,ny,nz,nmu,neta,nxi,nE);
Smat = get_mat_14D_fn(S,nx,ny,nz,nmu,neta,nxi,nE);
Fmat = get_mat_14D_fn(F,nx,ny,nz,nmu,neta,nxi,nE);
% Hmat = get_mat_10D_fn(H,nx,ny,nmu,neta,nE);
% Qvec = Q(:);
% Psivec = (Hmat-Smat-Fmat)\Qvec;
%% solve eigenvalue problem
A = (Hmat-Smat)\Fmat;
k = max(abs(eig(A)));

%% solve eigenvalue in tt
ttH = convert_tensor_op_to_tt_matrix_fn(H,14);
ttS = convert_tensor_op_to_tt_matrix_fn(S,14);
ttF = convert_tensor_op_to_tt_matrix_fn(F,14);
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
[ktt, ttPsi1] = tt_fixed_point_eig_solve(ttH, ttS, ttF, niter,...
  fixed_point_tol,tt_tol);

fprintf('TT eigenvalue = %.10f \n',ktt)
fprintf('Error in eigenvalue = %.2e \n',abs(ktt-k))
fprintf('compression ratio of H = %.5e \n', compress_ratio_tt(ttH));
fprintf('compression ratio of S = %.5e \n', compress_ratio_tt(ttS));
fprintf('compression ratio of F = %.5e \n', compress_ratio_tt(ttF));








