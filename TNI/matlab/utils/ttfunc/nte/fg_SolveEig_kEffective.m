function [ k, psi ] = fg_SolveEig_kEffective ( H, SigS, SigF, iV, L0S, ZL,...
    alpha, maxiter, epsi, print)

psi = ones( length(H), 1 );
k   = 1.0;

tic
for i = 1:maxiter

    % Save previous angular flux, calculate scalar flux and fission source.
    % Scalar Flux is phi = L0*S*psi where L0 is the matrix defined
    % Fission source is given by the product of the fission matrix and phi.
    psi_old      = psi;
    phi_old      = L0S*psi_old;
    fiss_src_old = L0S*SigF*ZL*phi_old;

    psi = ( H + alpha.*iV ) \ ( SigS*psi + (1/k).*SigF*psi );
    % phi = L0S*psi;
    % fiss_src = L0S*SigF*ZL*phi;

    k = k*sum( SigF*psi )/ sum( SigF*psi_old );

    lambda = sum( fiss_src ) / sum( fiss_src_old );
    maxdiff_phi = max( abs( (phi-phi_old)./phi_old ) );

    if ( print == true )
        fprintf( 'time = %f iter = %i k-eff eigenvalue = %f lambda-1 = %e max flux err = %e \n',...
            toc, i, k, abs(lambda-1), maxdiff_phi );
    end
    if ( abs(lambda-1) < epsi && maxdiff_phi < epsi )
        break
    end

end

return