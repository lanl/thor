close all; clear; clc;

%%
% f = @(x,y) 1 + sin(pi*x).*cos(pi*y);
f = @(x,y) sin(1.+2*pi*(x+y)) + 1/2*(x.^10 + y.^10);
%%
a = 0;
b = 1;

% Itrue = integral2(f,a,b,a,b);
Itrue = 1/11;
%%
Ns = [5,10,20];
Err = zeros(size(Ns));
for i = 1: numel(Ns)
  n = Ns(i);
  %% get quadrature points and weights
  x = linspace(0,1,n);
  dx = x(2) - x(1);
  % trapezoial weights
  w = ones(size(x));
  w(1) = 1/2;
  w(end) = 1/2;
  w = w*dx;
  %% Grid
  [Xgrid,Ygrid] = ndgrid(x);
  F = f(Xgrid,Ygrid); %evaluate function on grids

  %% Comtraction
  I = w*F*w';

  %%
  Err(i) = abs(I-Itrue);
  if i==1
    fprintf('n = %d - err = %.5e \n', n, Err(i));
  else
    convrate = log(Err(i)/Err(i-1))/log(1/2);
    fprintf('n = %d - err = %.5e - convrate = %.5e \n', n, Err(i),convrate);
  end
end

