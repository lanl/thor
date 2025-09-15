function I=applyQuadrature(tt,quad,h,dir)
    %
    % with 3 arguments, the output is a volume integral
    % with 4 arguments, the output is a surface intagral with h=1
    %
    if (nargin==4)
        h(:) = 1;
    else
        dir = -1;
    end 
    % 
    d = 2;
    % 
    r  = tt.r;
    n  = tt.n;
    %
    m  = quad.n;
    %
    I.n   = n/m;
    if(dir>0)
        I.n(dir) = n(dir);
    end
    I = create_tt(I.n,r,d);
    %
    for i=1:2
        %
        if(i~=dir)
            %
            Gk = reshape(tt{i},r(i),m,n(i)/m,r(i+1));
            Iq = zeros(r(i),n(i)/m,r(i+1));
            %
            for j=1:m
                Iq = Iq + reshape(Gk(:,j,:,:),r(i),n(i)/m,r(i+1))*quad.w(j);
            end
            %
            I{i} = h(i)*Iq;
            %
        else
            I{i} = tt{i};
        end
        %
    end
    %
end