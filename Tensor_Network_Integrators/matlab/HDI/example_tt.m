close all; clear; clc;

%%
% Testcase 1
f = @(x) sin(1.+2*pi*(sum(x))) + 1/2*(sum(x.^10));

%domain
a = 0;
b = 1;

Itrue = @(d) d/22; %exact value

%%
d = 20; %number of dimension
fprintf('Number of dimension = %d \n',d);
%%
Ns = [100,200,400]; %number of grid point per dimension

Err = zeros(size(Ns)); %preallocate error vector

for i = 1: numel(Ns)
  
  n = Ns(i);

  %% get quadrature points and weights
  % Trapazoidal rule
  x = linspace(0,1,n); %grid point
  dx = x(2) - x(1); % spatial step

  % trapezoial weights for 1D
  w = ones(size(x));
  w(1) = 1/2;
  w(end) = 1/2;
  w = w*dx;
  
  %% Approximate Ftt
  tol = 1e-12;
  Ftt = amen_cross2(n*ones(d,1),f,tol,x,'verb',0);
  
  %% FG checking
  % [Xgrid,Ygrid] = ndgrid(x);
  % F = f([Xgrid(:),Ygrid(:)]'); %evaluate function on grids
  % F = reshape(F,n*ones(1,d));
  % fprintf('TT err = %.5e \n',check_tt_rel_error(F,Ftt));
  %% Contraction
  % I = w*F*w';
  G = core2cell(Ftt); % get the core for the TT
  I = 1;
  for idim = 1:d
    %contract 1 core
    temp = tensorprod(G{idim},w,2,2);
    %accumuate to I
    I = I*temp;
  end
  %%
  Err(i) = abs(I-Itrue(d));
  if i==1
    fprintf('n = %d - err = %.5e \n', n, Err(i));
  else
    convrate = log(Err(i)/Err(i-1))/log(1/2);
    fprintf('n = %d - err = %.5e - convrate = %.5e \n', n, Err(i),convrate);
  end
end

