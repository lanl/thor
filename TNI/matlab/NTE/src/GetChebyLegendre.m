function [ mu, eta, xi, wgt ] = GetChebyLegendre ( N, idimen )

    [ x_tmp, w_tmp ] = ChebyshevLegendre( -1, 1, N );
    w_tmp = w_tmp / sum( w_tmp );

    kn = 0;
    nn = N/2;

    for ml = 1:nn
    
        z  = x_tmp(ml);
        wt = w_tmp(ml);

        if ( idimen == 2 )
            wt = 0.5*wt;
        elseif ( idimen == 3)
            wt = 0.25*wt;
        end

        sr1mz2 = sqrt( 1 - z^2 );

        mup = nn;

        for mp = 1:mup
            kn = kn + 1;
            eta(kn) = sr1mz2 * sin( ( 0.25*pi)*(2*mp-1) / mup );
            mu(kn)  = sr1mz2 * cos( ( 0.25*pi)*(2*mp-1) / mup );
            wgt(kn)   = wt/mup;
        end

    end

    xi = sqrt( 1 - mu.^2 - eta.^2 );

return