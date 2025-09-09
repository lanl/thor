close all; clear; clc;
run setup2DSlab.m

%% parameters 1D
% XY Dimensions
x0 = 0.0; x1 = 4.2;

ibl = 0;
ibr = 0;

idimen = 1;
octs   = 2^idimen;

% Load Cross Section Data
%load( 'EightGrpData.mat' );
load( 'OneGrpData.mat' );
% load( 'TwoGrpData.mat' );
G = length( sigt );

% Problem Parameters
M = 2;
L = 2;
J = 0;
K = 0;

ords = L;

% Fission Neutron Group PDF
chi    = ones( G,1 )/G;
chi    = chi.*ones( G,G );

% Quadrature Information
[ mu, wgt ] = AngularQuad1DSlab ( L, -1, 1 );

%Source
Q  = ones( G*(M+1)*ords, 1 );

%Include Fission Source
nosigf = 0;
nosigs = 0;

%%
% construct matrix from H
[ H, psi,Qbc] = get_fg_1D_operator_fn (sigt, sigs, nusigf, chi, vel,...
    mu, wgt, Q, x0, x1, M, ibl, ibr, nosigf, nosigs );

% [ psi_e, k, alpha ] = Solve1DSlab ( sigt, sigs, nusigf, chi, vel,...
%     mu, wgt, Q, x0, x1, M, ibl, ibr, nosigf, nosigs );

% construct matrix from Htt
nx = M+1;
nmu = L;
nE = G;
[ ttH] = get_tt_1D_operator_fn (sigt, x0, x1, L, M+1);
Q = 1*ones(nx,nmu,nE); %this is a vector of ones
Q(1,1,:) = 0; %BC mu pos
Q(nx,nmu,:) = 0; % BC mu neg
ttpsi = full(ttH)\Q(:);
% construct matrix from Hqtt