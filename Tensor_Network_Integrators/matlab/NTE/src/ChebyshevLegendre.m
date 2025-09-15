function [ x, w ] = ChebyshevLegendre ( a, b, n )

    eps = 3e-15;
    maxit = 10;

    m  = floor( (n + 1)/2 );
    xm = 0.5*( b + a );
    xl = 0.5*( b - a );

    num = 1:m;
    z = cos( pi*(num-0.25)/(n+0.5) )';

    x  = zeros( n, 1 );
    p1 = ones( m, 1 );
    p2 = zeros( m, 1 );
    p3 = zeros( m, 1 );
    pp = zeros( m, 1 );
    z1 = zeros( m, 1 );

    unfinished = ones( m, 1 );

    for its = 1:maxit
        p1( unfinished == 1 ) = 1;
        p2( unfinished == 1 ) = 0;
        for j = 1:n
            p3( unfinished == 1 ) = p2( unfinished == 1 );
            p2( unfinished == 1 ) = p1( unfinished == 1 );
            p1( unfinished == 1 ) = ((2*j-1).*z(unfinished==1).*p2( unfinished == 1 ) -...
                (j-1).*p3( unfinished == 1 ))./j;
        end
        for j = 1:length( unfinished )
            if unfinished(j) == 0
               continue
            end
            pp(j) = n * ( z(j)*p1(j)-p2(j) ) / (z(j)^2-1);
            z1(j) = z(j);
            z(j)  = z1(j) - p1(j) / pp(j);
            unfinished(j) = ( abs( z(j) - z1(j) ) > eps );
        end
        if ( ~any( unfinished ) )
            break
        end
    end

    x(1:m)        = xm - xl.*z;
    x(n:-1:n-m+1) = xm + xl.*z;
    w(1:m)        = 2.*xl ./ ( (1-z.^2).*pp.^2 );
    for its = n:-1:n-m+1
        w(its) = w( n-its+1 );
    end

return