function data=Q0(data)
    %
    gam = data.gam;
    %
    rho = zeros(size(data.X));
    v   = zeros(size(data.X));
    p   = zeros(size(data.X));
    %
    idxL = data.Y <0.5;
    idxR = data.Y>=0.5;
    %
    rho(idxL) = 2;
    p(idxL)   = 1+2*data.Y(idxL);
    %
    rho(idxR) = 1;
    p(idxR)   = 1.5+data.Y(idxR);
    %
    v(:) = -0.025*sqrt(gam*p(:)./rho(:)).*cos(8*pi*data.X(:));
    %
    %
    data.Q(1,:,:,:) = rho;
    data.Q(2,:,:,:) = 0;
    data.Q(3,:,:,:) = rho.*v;
    data.Q(4,:,:,:) = 0;
    data.Q(5,:,:,:) = p/(gam-1) + 0.5*rho.*(v.^2);
    %
end