function [R,L]=eigVecs(Q,d,gam)
    %
    Neq = size(Q,1);
    Nx  = size(Q,2);
    Ny  = size(Q,3);
    Nz  = size(Q,4);
    %
    R = zeros(Neq,Neq,Nx,Ny,Nz);
    L = zeros(Neq,Neq,Nx,Ny,Nz);
    %
    rho = Q(1,:,:,:);
    u   = Q(2,:,:,:)./rho;
    v   = Q(3,:,:,:)./rho;
    w   = Q(4,:,:,:)./rho;
    %
    ek = 0.5*(u.^2 + v.^2 + w.^2);
    %
    p = (gam-1)*(Q(5,:,:,:) - rho.*ek); 
    %
    c2 = gam*p./rho;
    c  = sqrt(c2);
    %
    nx = (1==d) + 0;
    ny = (2==d) + 0;
    nz = (3==d) + 0;
    %
    if (d==1)
        Vn=u;
    elseif (d==2)
        Vn=v;
    else
        Vn=w;
    end
    %
    phi = (gam-1)*ek;
    a1  = (gam-1);
    a2  = 1./(rho.*c*sqrt(2));
    a3  = rho./(c*sqrt(2));
    a4  = (phi + c2)/(gam-1);
    a5  = 1 - phi./c2;
    a6  = phi/(gam-1);
    %
    R(1,1,:,:,:) = nx;
    R(2,1,:,:,:) = nx*u;
    R(3,1,:,:,:) = nx*v + nz*rho;
    R(4,1,:,:,:) = nx*w - ny*rho;
    R(5,1,:,:,:) = nx*a6 + rho.*(nz*v - ny*w);
    %
    R(1,2,:,:,:) = ny;
    R(2,2,:,:,:) = ny*u - nz*rho;
    R(3,2,:,:,:) = ny*v;
    R(4,2,:,:,:) = ny*w + nx*rho;
    R(5,2,:,:,:) = ny*a6 + rho.*(nx*w - nz*u);
    %
    R(1,3,:,:,:) = nz;
    R(2,3,:,:,:) = nz*u + ny*rho;
    R(3,3,:,:,:) = nz*v - nx*rho;
    R(4,3,:,:,:) = nz*w;
    R(5,3,:,:,:) = nz*a6 + rho.*(ny*u - nx*v);
    %
    R(1,4,:,:,:) = a3;
    R(2,4,:,:,:) = a3.*(u + nx*c);
    R(3,4,:,:,:) = a3.*(v + ny*c);
    R(4,4,:,:,:) = a3.*(w + nz*c);
    R(5,4,:,:,:) = a3.*(a4 + c.*Vn);
    %
    R(1,5,:,:,:) = a3;
    R(2,5,:,:,:) = a3.*(u - nx*c);
    R(3,5,:,:,:) = a3.*(v - ny*c);
    R(4,5,:,:,:) = a3.*(w - nz*c);
    R(5,5,:,:,:) = a3.*(a4 - c.*Vn);
    %
    L(1,1,:,:,:) = nx*a5 - (nz*v-ny*w)./rho;
    L(2,1,:,:,:) = ny*a5 - (nx*w - nz*u)./rho;
    L(3,1,:,:,:) = nz*a5 - (ny*u - nx*v)./rho;
    L(4,1,:,:,:) =  a2.*(phi  - c.*Vn);
    L(5,1,:,:,:) =  a2.*(phi  + c.*Vn);
    %
    L(1,2,:,:,:) = nx*a1*u./c2;
    L(2,2,:,:,:) = ny*a1*u./c2 - nz./rho;
    L(3,2,:,:,:) = nz*a1*u./c2 + ny./rho;
    L(4,2,:,:,:) = -a2.*(a1*u - nx*c);
    L(5,2,:,:,:) = -a2.*(a1*u + nx*c);
    %
    L(1,3,:,:,:) = nx*a1*v./c2 + nz./rho;
    L(2,3,:,:,:) = ny*a1*v./c2;
    L(3,3,:,:,:) = nz*a1*v./c2 - nx./rho;
    L(4,3,:,:,:) = -a2.*(a1*v - ny*c);
    L(5,3,:,:,:) = -a2.*(a1*v + ny*c);
    %
    L(1,4,:,:,:) = nx*a1*w./c2 - ny./rho;
    L(2,4,:,:,:) = ny*a1*w./c2 + nx./rho;
    L(3,4,:,:,:) = nz*a1*w./c2;
    L(4,4,:,:,:) = -a2.*(a1*w - nz*c);
    L(5,4,:,:,:) = -a2.*(a1*w + nz*c);
    %
    L(1,5,:,:,:) = -nx*a1./c2;
    L(2,5,:,:,:) = -ny*a1./c2;
    L(3,5,:,:,:) = -nz*a1./c2;
    L(4,5,:,:,:) =  a1*a2;
    L(5,5,:,:,:) =  a1*a2;
    %
%   Eig = zeros(size(R));
%   %%
%   Eig(1,1,:,:,:) = Vn; 
%   Eig(2,2,:,:,:) = Vn; 
%   Eig(3,3,:,:,:) = Vn; 
%   Eig(4,4,:,:,:) = Vn + c; 
%   Eig(5,5,:,:,:) = Vn - c; 
    %
%    J = zeros(size(L));
%    %
%    b1 = gam*(Q(5,:,:,:)./Q(1,:,:,:)) - phi;
%    b2 = gam-1;
%    b3 = gam-2;
%    %
%    %J(1,1,:,:,:) = 0;
%    J(1,2,:,:,:) = nx;
%    J(1,3,:,:,:) = ny;
%    J(1,4,:,:,:) = nz;
%    %J(1,5,:,:,:) = 0;
%    %
%    J(2,1,:,:,:) = nx*phi - u.*Vn;
%    J(2,2,:,:,:) = Vn - b3*nx*u;
%    J(2,3,:,:,:) = ny*u - b2*nx*v;
%    J(2,4,:,:,:) = nz*u - b2*nx*w;
%    J(2,5,:,:,:) = b2*nx;
%    %
%    J(3,1,:,:,:) = ny*phi - v.*Vn;
%    J(3,2,:,:,:) = nx*v - b2*ny*u;
%    J(3,3,:,:,:) = Vn - b3*ny*v;
%    J(3,4,:,:,:) = nz*v - b2*ny*w;
%    J(3,5,:,:,:) = b2*ny;
%    %
%    J(4,1,:,:,:) = nz*phi - w.*Vn;
%    J(4,2,:,:,:) = nx*w - b2*nz*u;
%    J(4,3,:,:,:) = ny*w - b2*nz*v;
%    J(4,4,:,:,:) = Vn - b3*nz*w;
%    J(4,5,:,:,:) = b2*nz;
%    %
%    J(5,1,:,:,:) = Vn.*(phi-b1);
%    J(5,2,:,:,:) = nx*b1 - b2*u.*Vn;
%    J(5,3,:,:,:) = ny*b1 - b2*v.*Vn;
%    J(5,4,:,:,:) = nz*b1 - b2*w.*Vn;
%    J(5,5,:,:,:) = gam*Vn;