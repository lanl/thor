function I=applyQuadrature(F,data,t)
    %
    Neq = data.Neq;
    Nx  = data.Nx;
    Ny  = data.Ny;
    %
    vol = data.vol;
    %
    wij = data.quad.wij;
    rx  = data.quad.rx;
    ry  = data.quad.ry;
    %
    I = zeros(Neq,Nx,Ny);
    %
    for j = 1:Ny
        for i=1:Nx    
            %
            qx = (rx + i - 1)*data.dx;
            qy = (ry + j - 1)*data.dy;
            %
            if (nargin==3)
                Int = F(qx,qy,t,data);
            else
                Int = F(qx,qy,data);
            end
            %
            for eq=1:Neq
                I(eq,i,j) = sum(wij(:)'.*Int(eq,:)*vol,"all");
            end
        end
    end
    %
end