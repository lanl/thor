close all; clear; clc;

% close all; clear; clc;
%main file to test Gauss Quadrature
run setup.m
%% define function
if 0
  d=2;
  f1 = @(x) (x-1/2).^10.*sign(x-1/2);
  T10 = @(x) 512.*x.^10 - 1280.*x.^8 + 1120.*x.^6 - 400.*x.^4 + 50.*x.^2 - 1;
  
  f = @(x,y) f1(x) + (chebyshevT(10,x) + chebyshevT(10,y))/d;
  f = @(x,y) f1(x) + (T10(x) +T10(y))/d;
  Itrue = -1/99;
end

%%
[f,Itrue] = select_integrand_function(7);

%% compute quadrature

a = 0;
b = 1;
ncell = 2;
dx = (b-a)/ncell;
Ns = [1:10];

for norder = 1:10
  
  xgrid = zeros(1,norder*ncell);
  w = zeros(size(xgrid));
  
  for i = 1:ncell
    tempx0 = a + (i-1)*dx;
    tempx1 = tempx0 + dx;
    
    [tempx,tempw] = legpts(norder,[tempx0,tempx1]); %get GL nodes and weights
    xgrid((i-1)*norder+1:i*norder) = tempx;
    w((i-1)*norder+1:i*norder) = tempw;
  end
  
  
  %evaluate function
  
  [X,Y] = meshgrid(xgrid,xgrid);
  feval = f(X,Y);
  
  %compute integral
  
  I = w*feval*w';
  fprintf('norder = %d, I = %.5e, Err = %.5e\n', norder, I, abs(I-Itrue));
  % Err = mp(abs(w*feval*w'-Itrue));
end
