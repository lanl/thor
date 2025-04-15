function [result] = tt_GL_fn(itest, d, ncell, method, norder, tol)
%% define function
[f,Itrue] = select_integrand_function(itest,d);

%% compute quadrature
total_time_start = tic;

a = 0;
b = 1;
% ncell = 10;
dx = (b-a)/ncell;
% norder = 4;

xgrid = zeros(1,norder*ncell);
w = zeros(size(xgrid));

% get weight and xgrid
for i = 1:ncell
  
  tempx0 = a + (i-1)*dx;
  tempx1 = tempx0 + dx;
  
  %% get node and weight using chebfun
  switch method
    case 'GL'
      [tempx,tempw] = legpts(norder,[tempx0,tempx1]);
    case 'CC'
      [tempx,tempw] = chebpts(norder,[tempx0,tempx1]);
    case 'Simpson'
      wq = [  1/3, 4/3, 1/3 ];
      sq = [-1.0,   0.0,   1.0];
      tempx = tempx0 + (tempx1-tempx0)*(1+sq)/2;
      tempw = (tempx1-tempx0)*wq/2;
    case 'Trapezoidal'
      wq = [1,1];
      sq = [-1,1];
      tempx = tempx0 + (tempx1-tempx0)*(1+sq)/2;
      tempw = (tempx1-tempx0)*wq/2;
  end
  %   [tempx,tempw] = gauss_legendre_quadrature(norder,tempx0,tempx1);
  %%
  xgrid((i-1)*norder+1:i*norder) = tempx;
  w((i-1)*norder+1:i*norder) = tempw;
end

%% evaluate function
n = numel(xgrid);
Ns = ones(1,d)*n;


time_cross_start = tic;
ftt = amen_cross2(Ns, f, tol, xgrid);
result.time_cross = toc(time_cross_start);

%% dot product with W to compute the quadrature
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
result.time = toc(total_time_start);
result.r = ftt.r;
result.n = ftt.n;
result.Err = abs(S-Itrue);
result.tol = tol;

fprintf('Err = %.5e, time = %.5e \n', result.Err, result.time);
end


