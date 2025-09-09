function Err = full_GL_2D_fn(ncell)

% close all; clear; clc;
%main file to test Gauss Quadrature

%% define function
q = 30;
d = 2;
f= @(x) mp((q+1)/d.*(x.^q + y.^q));
Itrue = 1;
%% compute quadrature
q = mp(10);
a = mp(0);
b = mp(1);
% ncell = 10;
dx = mp((b-a)/ncell);
norder = 4;

xgrid = mp(zeros(1,norder*ncell));
w = mp(zeros(size(xgrid)));

for i = 1:ncell
  tempx0 = mp(a + (i-1)*dx);
  tempx1 = mp(tempx0 + dx);
  
  [tempx,tempw] = mp.GaussLegendre(norder,tempx0,tempx1);
  xgrid((i-1)*norder+1:i*norder) = tempx;
  w((i-1)*norder+1:i*norder) = tempw;
end


%evaluate function
[X,Y] = meshgrid(xgrid,xgrid);
feval = mp(f(X,Y,q));

Err = mp(abs(w*feval*w'-Itrue));
fprintf('Err = %.50e \n',Err);


end


