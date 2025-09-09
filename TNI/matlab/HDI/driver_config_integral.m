%% Compute the configuration integrals
close all; clear; clc;
run setup.m;

%%
tol=1e-2;

%% setup the function
%parameters
a = 1;
param_d = 1;
param_e = 1.2; %temperature --- range from very small to 2
x0 = 0;
x1 = 1;
N = 20;  % number of particles 
d = 3*N; %number of dimension 
scale_factor = 1;
f = @(x) scale_factor*(Zn(x,N, param_d, param_e)); %without d^3N factor
% f = @g2;
%% minimize u
% U = @(x) -u_func(x(1:3),x(4:6),param_d,param_e);
% x0 = rand(1,6);
% lb = zeros(1,6);
% ub = ones(1,6);
% [xmin, fmin] = fmincon(U, x0,[],[],[],[],lb,ub);

%% get weights and nodes
ncell = 10;
norder = 2;

dx = (x1-x0)/ncell;

xgrid = zeros(1,norder*ncell);
w = zeros(size(xgrid));

for i = 1:ncell
  tempx0 = x0 + (i-1)*dx;
  tempx1 = tempx0 + dx;
  %get GL weight and node
  [tempx,tempw] = legpts(norder,[tempx0,tempx1]);
  xgrid((i-1)*norder+1:i*norder) = tempx;
  w((i-1)*norder+1:i*norder) = tempw;
end

%% tt evaluate the function
tic
n = numel(xgrid);
Ns = ones(1,d)*n;
ftt = amen_cross2(Ns, f, tol, xgrid);
toc
%% dot product to compute the integrals
G = core2cell(ftt);
S = 1;
for i = 1:d
  if i~=d
    g = reshape(permute(G{i},[2,1,3]),n,[]);
    S = S*reshape((w*g),size(G{i},1),[]);
    %reshape to matrix
  else %last core
    S = S*(G{d}*w');
  end
end


fprintf('N = %d particles, a = %.5f, d = %.5f, e = %.5f \n',N, a, param_d, param_e);
fprintf('Zn = %.5f \n', S + 1);

%% subroutines
function f = Zn(x, N, d, e)

%compute g
g = 0;
for p1 = 1:(N-1)
  for p2 = (p1+1) : N
    g = g + u_func(x(3*(p1-1)+1:3*p1), x(3*(p2-1)+1:3*p2), d, e);
  end
end
f = exp(-g);

end

function u = u_func(x, y, d, e)
% compute the u(rij) function
r = norm(x-y);
if r == 0
  r = 1e-16;
end
u = 4*e*((d/r)^12 - (d/r)^6);
end
