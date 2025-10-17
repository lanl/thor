function N0n = GetDNUBSf (knot,Zita,m)
KnotSpan  = GetSpan(knot,Zita);
[n,p]     = EvaluateKnot(knot);
N0n       = zeros(p+1,m+1);
iKnotSpan = KnotSpan;
iXi       = Zita;
Ni        = zeros(p+1,m+1);
ndu       = zeros(p+1,p+1);
left      = zeros(p+1);
right     = zeros(p+1);
a         = zeros(2,p+1);
ndu(1, 1) = 1;


    for j = 1:p
        left(j+1)  = iXi - knot(iKnotSpan+1-j);
        right(j+1) = knot(iKnotSpan+j) - iXi;
        saved = 0;
        for r = 0:j-1
            % Lower triangle
            ndu(j+1,r+1) = right(r+2)+left(j-r+1);
            temp         = ndu(r+1,j)/ndu(j+1,r+1);
            % Upper triangle
            ndu(r+1,j+1) = saved+right(r+2)*temp;
            saved = left(j-r+1)*temp;
        end
        ndu(j+1,j+1) = saved;
    end
    % Load the basis function
    Ni(:,1) = ndu(:,p+1);
    % Compute the derivatives
    for r =0:p % Loop over function index
        s1 = 0;
        s2 = 1; % Alternate rows in array a
        a(1,1) = 1;
        % Loop to compute kth derivative
        for k =1:m
            d  =0;
            rk =r-k;
            pk =p-k;
            if (r >= k)
                a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1);
                d         = a(s2+1,1)*ndu(rk+1,pk+1);
            end
            if (rk>=-1)
                j1=1;
            else
                j1=-rk;
            end
            if ((r-1) <= pk)
                j2=k-1;
            else
                j2=p-r;
            end
            for j=j1:j2
                a(s2+1,j+1)=(a(s1+1,j+1)-a(s1+1,j))/ndu(pk+2,rk+j+1);
                         d =d+a(s2+1,j+1)*ndu(rk+j+1,pk+1);
            end
            if (r<= pk)
                a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
                d           = d+a(s2+1,k+1)*ndu(r+1,pk+1);
            end
            Ni(r+1,k+1) = d;
            j=s1;
            s1=s2;
            s2=j;
        end
    end
    r = p;
    for k = 1:m
        Ni(:,k+1) = Ni(:,k+1)*r;
        r = r*(p-k);
    end
    N0n(:, :) = Ni;
end