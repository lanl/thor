clc, clear

% Slab Dimensions
a = 0.0;
b = 1.0;

% Total Cross Section
sigt = 0.45391;
% sigs = 0.34612;
sigs = 0.1184;
nusigf = 0.24683;
mu   = [ -0.57735, 0.57735 ];
wgt  = [ 0.5, 0.5 ];

sigt   = [ 2.3336e-01 2.3105e-01 3.0084e-01 2.8809e-01 3.8538e-01 5.0120e-01 6.1661e-01 7.8866e-01];
sigs   = load( './Data/sigma_s.txt' );
nusigf = load( './Data/eig_nusigma_f.txt' );
chi    = load( './Data/eig_chi.txt' );
nusigf = chi.*nusigf;

% Problem Parameters
M = 2;
L = 2;
G = 8;

% Generate identity matrix for L angles
IL = eye( L );
IM = eye( M );
IG = eye( G );

% Generate Projection Matrices ( edges <-> cell-center )
Z = zeros(M,M);
Zm = Z(M,:)';
Z = [eye(M,M); zeros(1,M)];
Zb = zeros(M+1,1); Zb(end) = 1;

% Generate Interaction Term
Sm = 1.0.*( [Zm,diag( ones( M,1 ),0 ) ] );
Sp = 1.0.*( [diag( ones( M,1 ),0),Zm ] );
Sm = kron( IG, kron( IL, 0.5*Sm ) );
Sp = kron( IG, kron( IL, 0.5*Sp ) );
S = Sp + Sm;

% Generate Total Cross Section Matrix
Sigma = zeros( G*L*M, G*L*M );
for g = 1:G
    for l = 1:L
        for m = 1:M
            Sigma = Sigma + sigt(g) * kron( IG(:,g), kron( IL(:,l), IM(:,m) ) ) *...
                kron( IG(:,g), kron( IL(:,l), IM(:,m) ) )';
        end
    end
end

% Generate Scattering Cross Section Matrix
SigmaS = zeros( G*M, G*M );
for gp = 1:G
    for g = 1:G
        for m = 1:M
            SigmaS = SigmaS + sigs(g,gp) * kron( IG(:,g), IM(:,m) ) *...
                kron( IG(:,gp), IM(:,m) )';
        end
    end
end

% Generate Fission Cross Section Matrix
%SigmaF = diag( nusigf.*ones( M,1 ) );
SigmaF = zeros( G*L*M, G*L*M );
for g = 1:G
  for l = 1:L
    for m = 1:M
      SigmaF = SigmaF + nusigf(g) * kron( IG(:,g), kron( IL(:,l), IM(:,m) ) ) *...
        kron( IG(:,g), kron( IL(:,l), IM(:,m) ) )';
    end
  end
end

L0  = kron( IG, kron( [1 1]*diag( wgt ), IM ) );
L0p = kron( IG, kron( ones(L,1), IM ) );

% Generate Derivative Matrix
dx = ( b -a )/M;
dx = diag( dx*ones( 1,M ) );
Dm = -1.0.*( [diag( ones( M,1 ),0),Zm ] );
Dp = 1.0.*( [Zm,diag( ones( M,1 ),0 ) ] );
D  = Dm + Dp;

% Generate Projection Matrices ( edges <-> cell-center ) for L angles
Z  = kron( IG, kron( IL, Z ) );
ZB = kron( IG, kron( IL, Zb ) );

% Generate Boundary Conditon Matrix
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
B = kron( IG, B );

% Generate matrix H with size L(M+1) x L(M+1)
Hmu  = Z * ( kron( IG, kron( diag(mu), dx\D ) ) );
Hmup = Z * ( kron( IG, kron( diag(mu), dx\Dp ) ) );
Hmum = Z * ( kron( IG, kron( diag(mu), dx\Dm ) ) );

Hsig = Z * Sigma * S;
Hsigm = Z * Sigma * Sm;
Hsigp = Z * Sigma * Sp;

SigmaS = Z*L0p*SigmaS*L0*S;
%SigmaF = Z*L0p*SigmaF*L0*S;
SigmaF = Z*SigmaF*S;

Hbc  = ZB*B;

H = Hmu + Hsig + Hbc;

%SigS = Z*L0p*SigmaS*L0*S;
%SigF = Z*L0p*SigmaF*L0*S;

% External source ones, apply vacuum BCs
Q = ones( G*(M+1)*L,1 );
Qbc = Q - ZB*B*Q;

% Generate Scattering Matrix

% MATLAB inversion of H
psi = ( H - SigmaF ) \ ( Qbc );
%psi = (H-SigS-SigF) \ ( Qbc );

% Matrix Splitting
%A = Hmup + Hmum + Hsigm + Hbc;
%B = -Hsigp;
A = Hmum + Hmup + Hbc;
B = -Hsigp - Hsigm + SigmaF;
%B = -Hsigp - Hsigm + SigS + SigF;


% Fixed Point Solve, Need to Add Error Check and Test More Splittings
e = zeros(100);
psi0 = ones( G*L*(M+1), 1 );
for i = 1:1000
   psi_old = psi0;
   psi0 = A\B*psi0 + A\Qbc;
   e(i) = norm( psi0 - psi_old );
end