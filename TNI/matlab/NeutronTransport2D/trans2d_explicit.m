close all; clear; clc;
run setup2DSlab.m

%% Slab Dimensions
x0 = 0.0; x1 = 4.2;
y0 = 0.0; y1 = 4.2;

% Total Cross Section
sigt = 0.45391;
sigt = 0;
% sigt   = [ 2.3336e-01 2.3105e-01 3.0084e-01 2.8809e-01 3.8538e-01 5.0120e-01 6.1661e-01 7.8866e-01];
%sigt   = [ 2.9138e-1 4.7284e-1 ];
sigs = 0.34612;
nusigf = 0.24683;

mu   = [ -0.57735027, 0.57735027, -0.57735027, 0.57735027 ];
eta  = [ -0.57735027, -0.57735027, 0.57735027, 0.57735027 ];

wgt  = [ 0.25, 0.25, 0.25, 0.25 ];
% mu = 0*mu;
% eta = 0*eta;
% Problem Parameters
M = 2; %x
J = 2; %y
L = 2; % mu
N = 2; % eta
G = length( sigt );

chi    = ones( G,1 )/G;
nusigf = chi.*nusigf;


%% Generate identity matrix for L angles
IL = eye( L );
IM = eye( M );
IG = eye( G );
ILN = eye( L^2 ); %for double angle

%% Generate Projection Matrices ( edges <-> cell-center )
ZM = zeros(M,M); ZN = zeros(J,J);
ZMm = ZM(M,:)'; ZNm = ZN(J,:)';
ZM = [eye(M,M); zeros(1,M)]; ZN = [eye(J,J); zeros(1,J)];
ZMb = zeros(M+1,1); ZMb(end) = 1; ZNb = zeros(J+1,1); ZNb(end) = 1;

Z = [eye(M*J); zeros((M+1)*(J+1)-M*J,M*J)];
II = eye( (M+1)*(J+1)-M*J );
Zb = [ zeros( M*J, (M+1)*(J+1)-M*J); II ];

%% Generate Interaction Term
SM = 0.5.*( [ZMm,diag( ones( M,1 ),0 ) ] ) + 0.5.*( [diag( ones( M,1 ),0),ZMm ] );
SJ = 0.5.*( [ZNm,diag( ones( J,1 ),0 ) ] ) + 0.5.*( [diag( ones( J,1 ),0),ZNm ] );
S = kron( SJ, SM );

% Generate Derivative Matrix d/dx
dx = ( x1 - x0 )/M;
dx = diag( dx*ones( 1,M ) );
DM = -1.0.*( [diag( ones( M,1 ),0),ZMm ] ) + 1.0.*( [ZMm,diag( ones( M,1 ),0 ) ] );

% Generate Derivative Matrix d/dy
dy = ( y1 - y0 )/J;
dy = diag( dy*ones( 1,J ) );
DN = -1.0.*( [diag( ones( J,1 ),0),ZNm ] ) + 1.0.*( [ZNm,diag( ones( J,1 ),0 ) ] );

%% Sigma = kron( eye(L*N), sigt.*eye(M*J) );
for g = 1:G
  Sigma{g} = kron( ILN, sigt(g).*eye(M*J) );
end
Sigma = blkdiag(Sigma{:});

% Generate Boundary Condition Matrices
em = eye(M+1);
e0m = em(1,:)'; emm = em(M+1,:)';
en = eye(J+1);
e0n = en(1,:)'; enn = en(J+1,:)';

% These are functions of the sign of mu and eta. They are in the following
% order:
% 1) mu<0, eta<0
% 2) mu>0, eta<0
% 3) mu<0, eta>0
% 4) mu>0, eta>0
E{1} = [ kron( enn', eye(M+1) ); kron( [ eye(J), zeros(J,1) ], emm' ) ];
E{2} = [ kron( enn', eye(M+1) ); kron( [ eye(J), zeros(J,1) ], e0m' ) ];
E{3} = [ kron( e0n', eye(M+1) ); kron( [ zeros(J,1), eye(J) ], emm' ) ];
E{4} = [ kron( e0n', eye(M+1) ); kron( [ zeros(J,1), eye(J) ], e0m' ) ];

Cx = kron( SJ, dx\DM );
Cy = kron( dy\DN, SM );

B  = kron( IG, blkdiag(E{:}) );
Z  = kron( IG, kron( ILN, Z ) );
S  = kron( IG, kron( ILN, S ) );
ZB = kron( IG, kron( ILN, Zb ) );
Cx = kron( IG, kron( diag(mu), Cx ) );
Cy = kron( IG, kron( diag(eta), Cy ) );
Q  = ones( G*(M+1)*(J+1)*(L^2), 1 );

size(ZB*B)

H   = Z*(Cx+Cy+Sigma*S) + ZB*B;
QQ  = Q - ZB*B*Q;

%% full-grid
psi_edge_centered = H\QQ;
psi_cell_centered = Z'*(Z*S + ZB*B)*psi_edge_centered;
%% tt
tol = 1e-6;
Hs = reshape(H,[G,L,N,(M+1),(J+1),G,L,N,(M+1),(J+1)]);
Htt = tt_matrix(Hs,tol);

%%
QQs = reshape(QQ,[G,L,N,(M+1),(J+1)]);
QQtt = tt_tensor(QQs);

norm(QQ - full(QQtt))

%%
psitt = amen_solve2(Htt,QQtt,1e-6,{'nswp',5});

fprintf('tt error in Psi = %.2e \n', check_tt_error(psi_edge_centered,psitt));



