close all; clear; clc;

% Slab Dimensions
a = 0.0;
b = 1.0;

% Total Cross Section
sigt = 0.45391;
mu   = [ -0.57735, 0.57735 ];

% Problem Parameters
M = 2;
L = 2;
G = 1;

% Generate identity matrix for L angles
IL = eye( L );

% Generate Projection Matrices ( edges <-> cell-center )
Z = zeros(M,M);
Zm = Z(M,:)';
Z = [eye(M,M); zeros(1,M)];
Zb = zeros(M+1,1); Zb(end) = 1;

%% Generate Interaction Term
Sm = 1.0.*( [Zm,diag( ones( M,1 ),0 ) ] );
Sp = 1.0.*( [diag( ones( M,1 ),0),Zm ] );
Sm = kron( IL, 0.5*Sm );
Sp = kron( IL, 0.5*Sp );
S = Sp + Sm;

%% Generate Total Cross Section Matrix
Sigma = kron( IL, diag( sigt.*ones( M,1 ) ) );

% Generate Derivative Matrix
dx = ( b -a )/M;
dx = diag( dx*ones( 1,M ) );
Dm = -1.0.*( [diag( ones( M,1 ),0),Zm ] );
Dp = 1.0.*( [Zm,diag( ones( M,1 ),0 ) ] );
D  = Dm + Dp;

%% Generate Projection Matrices ( edges <-> cell-center ) for L angles
Z = kron( IL, Z );
ZB = kron( IL, Zb );

%% Generate Boundary Conditon Matrix
e = eye(M+1);
e0m = e(1,:)';
emm = e(M+1,:)';
B1 = e0m';
B2 = emm';

for i = 1:length(mu)
    if ( mu(i) < 0 )
        B{i} = B2;
    else
        B{i} = B1;
    end
end

B = blkdiag(B{:});

%% Generate matrix H with size L(M+1) x L(M+1)
Hmu  = Z * ( kron( diag(mu), dx\D ) );
Hmup = Z * ( kron( diag(mu), dx\Dp ) );
Hmum = Z * ( kron( diag(mu), dx\Dm ) );

Hsig = Z * Sigma * S;
Hsigm = Z * Sigma * Sm;
Hsigp = Z * Sigma * Sp;

Hbc  = ZB*B;

H = Hmu + Hsig + Hbc;

% External source ones, apply vacuum BCs
Q = ones( (M+1)*L,1 );
Qbc = Q - ZB*B*Q;

% MATLAB inversion of H
psi = H \ ( Qbc );

% Matrix Splitting
A = Hmup + Hmum + Hsigm + Hbc;
B = -Hsigp;

% Fixed Point Solve, Need to Add Error Check and Test More Splittings
psi0 = ones( L*(M+1), 1 );
for i = 1:100
   psi_old = psi0;
   psi0 = A\B*psi0 + A\Qbc
end