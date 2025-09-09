function data=recon5(data,d)
    % get size
    Neq  = data.tt.Neq;
    %
    Qstencil_Neq = getWeno5StencilNeq(data.tt.Q,d);
    %
    % loop over equations
    %
    for i=1:Neq
        %
        Qvar = Qstencil_Neq{i+2*Neq};
        %
        if(norm(Qvar)>0)
            %
            Qstencil = cell(1,5);
            %
            for j=1:5
                Qstencil{j} = Qstencil_Neq{i+Neq*(j-1)};
            end
            %
            QR = cross_interpolation(Qstencil,  @(x) data.recon5fun(x,+1,data.weno.weps), data.tt.eps_cr, Qstencil{3}, 0);
            QL = cross_interpolation(Qstencil,  @(x) data.recon5fun(x,-1,data.weno.weps), data.tt.eps_cr, Qstencil{3}, 0);
            %
            QR = recon5GQ(QR,data,d);
            QL = recon5GQ(QL,data,d);
    
            QR = round(QR,data.tt.eps);
            QL = round(QL,data.tt.eps);
            %
        else
            %
            QR = quad_tt(Qvar.n,Qvar.r,data.quad.n,d);
            QL = quad_tt(Qvar.n,Qvar.r,data.quad.n,d);
            %
        end
        %
        data.tt.QR{i} = QR;
        data.tt.QL{i} = QL;
        %
    end
    %
end
%
function tt=quad_tt(N,r,nq,dir)
    %
    tt     = tt_tensor;
    %
    tt.n      = N*nq;
    tt.n(dir) = N(dir);
    tt.d      = 2;
    tt.r      = r;
    tt.ps     = cumsum([1;tt.n.*tt.r(1:tt.d).*tt.r(2:tt.d+1)]);
    tt.core   = zeros(tt.ps(tt.d+1)-1,1);
    %
end