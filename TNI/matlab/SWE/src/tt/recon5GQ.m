function Qq=recon5GQ(Qq,data,dir)
    %
    if(dir==1)
        d = 2;
    elseif(dir==2)
        d = 1;
    end
    %
    Qstencil = getWeno5Stencil(Qq,d);
    %
    Qrecon = cross_interpolation(Qstencil,  @(x) data.recon5GQfun(x,data.weno,data.quad.n), data.tt.eps_cr, Qstencil{3}, 0);
    %
    r      = Qrecon.r;
    r(end) = 1;
    %
    Qq = quad_tt(Qrecon.n,r,data.quad.n,d);
    % 
    for j=1:data.quad.n
        %
        Qquad = quad_tt(Qrecon.n,r,data.quad.n,d);
        %
        jdx = j:data.quad.n:Qquad.n(d);
        %
        if(d==2)
            %
            Qquad{1} = Qrecon{1};
            %
            core2 = Qquad{2};
            %
            core2(:,jdx,:) = reshape(Qrecon{2}(:,:,j),[r(2),Qrecon.n(2),1]);
            %
            Qquad{2} = core2;
            %
        else
            %
            core1 = Qquad{1}; 
            core2 = Qquad{2};
            %
            core1(:,jdx,:) = Qrecon{1};
            core2(:,:,:)   = reshape(Qrecon{2}(:,:,j),[r(2),Qrecon.n(2),1]);
            %
            Qquad{1} = core1;
            Qquad{2} = core2;
            %
        end
        %
        Qq = round(Qq + Qquad,data.tt.eps);
        %
    end
    %
    return
    %
end
%
function tt = quad_tt(N,r,nq,dir)
    %
    tt = tt_tensor;
    %
    tt.n      = N;
    tt.n(dir) = N(dir)*nq;
    tt.d      = 2;
    tt.r      = r;
    tt.ps     = cumsum([1;tt.n.*tt.r(1:tt.d).*tt.r(2:tt.d+1)]);
    tt.core   = zeros(tt.ps(tt.d+1)-1,1);
    %
end