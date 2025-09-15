function [ psi, k] = Solve1DSlab ( sigt, sigs, nusigf, chi, vel,...
    mu, wgt, Q, a, b, M, ibl, ibr, nosigf, nosigs )

L = length( mu );
G = length( sigt );

%nusigf = chi*diag(nusigf);

% Generate identity matrix for L angles
IL = speye( L );
IM = speye( M );
IG = speye( G );

% Generate Projection Matrices ( edges <-> cell-center )
Z = sparse(M,M);
Zm = Z(M,:)';
Z = [speye(M,M); sparse(1,M)];
Zb = sparse(M+1,1); Zb(end) = 1;

% Generate Interaction Term
Sm = 1.0.*( [ Zm, diag( ones( M,1 ), 0 ) ] );
Sp = 1.0.*( [ diag( ones( M,1 ), 0 ), Zm ] );
S = 0.5*Sm + 0.5*Sp;

% Generate Total Cross Section Matrix
Sigma = sparse( G*L*(M+1), G*L*(M+1) );
for g = 1:G
    for m = 1:M
        Sigma = Sigma + sigt(g) * ...
            ( kron( IG(:,g)*IG(:,g)',...
              kron( IL, Z*(IM(:,m)*IM(:,m)')*S ) ) );
    end
end

% Generate Inverse Velocity Matrix
iVel = sparse( G*L*(M+1), G*L*(M+1) );
for g = 1:G
    for m = 1:M
        iVel = iVel + (1/vel(g)) * ...
            ( kron( IG(:,g)*IG(:,g)',...
              kron( IL, Z*(IM(:,m)*IM(:,m)')*S ) ) );
    end
end

% Generate Scattering Cross Section Matrix
SigmaS = sparse( G*M, G*M );
for gp = 1:G
    for g = 1:G
        for m = 1:M
            SigmaS = SigmaS + sigs(g,gp) * ...
                ( kron( IG(:,g)*IG(:,gp)', IM(:,m)*IM(:,m)' ) );
        end
    end
end

if ( nosigs == 1 )
    SigmaS(:) = 0;
end

% Generate Fission Cross Section Matrix
SigmaF = sparse( G*M, G*M );
for gp = 1:G
    for g = 1:G
        for m = 1:M
            SigmaF = SigmaF + chi(g,gp)*nusigf(gp) *...
                kron( IG(:,g)*IG(:,gp)', IM(:,m)*IM(:,m)' );
        end
    end
end

if ( nosigf == 1 )
    SigmaF(:) = 0;
end

% Generate Derivative Matrix
dx = ( b - a )/M;
dx = diag( dx*ones( 1,M ) );
Dm = -1.0.*( [ diag( ones( M,1 ), 0 ), Zm ] );
Dp = 1.0.*( [ Zm, diag( ones( M,1 ), 0 ) ] );
D  = Dm + Dp;

% Generate Projection Matrices ( edges <-> cell-center ) for L angles
ZB = kron( IG, kron( IL, Zb ) );

% Generate Boundary Conditon Matrix (Vacuum and Reflective)
e   = eye(M+1);
e0m = e(1,:)';
emm = e(M+1,:)';
B1  = e0m';
B2  = emm';

BL0 = sparse( G*L*(M+1), G*L*(M+1) );
BL1 = sparse( G*L*(M+1), G*L*(M+1) );
BR0 = sparse( G*L*(M+1), G*L*(M+1) );
BR1 = sparse( G*L*(M+1), G*L*(M+1) );

for l = 1:L/2
    BL0 = BL0 + kron( IG, kron( IL(:,l)*IL(:,l)', Zb*B2 ) );
    BR0 = BR0 + ibr.*kron( IG, kron( IL(:,l)*IL(:,L+1-l)', Zb*B2 ) );
end

for l=(L/2)+1:L
    BL1 = BL1 + kron( IG, kron( IL(:,l)*IL(:,l)', Zb*B1 ) );
    BR1 = BR1 + ibl.*kron( IG, kron( IL(:,l)*IL(:,L+1-l)', Zb*B1 ) );
end
BCV = BL0 + BL1;
BCR = BR0 + BR1;

% Generate matrix H with size L(M+1) x L(M+1)
Hmu = kron( IG, kron( diag(mu), Z*(dx\D) ) );

ZL0p = kron( IG, kron( ones( L, 1 ), Z ) );
L0S  = kron( IG, kron( ones( 1, L )*diag( wgt ), S ) );

SigS = ZL0p*SigmaS*L0S;
SigF = ZL0p*SigmaF*L0S;

H  = Hmu + Sigma + BCV - BCR;
iV = iVel + BCV - BCR;

% External source ones, apply vacuum BCs
Qbc = Q - BCV*Q;

% Solve Source Problem
psi = ( H - SigS - SigF ) \ ( Qbc );
psi = psi/norm(psi);
% k-effective Eigenvalue for System
A = ( H - SigS ) \ SigF;
k = max( abs( eigs( A ) ) );

% Alpha-eigenvalue for System
% A = SigS + SigF - H;
% B = iV;
% D = eigs( A, B );
% alpha = max( real( D( ~isinf( D ) ) ) );

return