function Err = tt_GL_fn_poly(d,q,norder, ncell, tol)
%tt cross without mp

%% define function
f= @(x) (q+1)/d.*sum(x.^q);
Itrue = 1;


%% compute quadrature
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
  
  [tempx,tempw] = gauss_legendre_quadrature(norder,tempx0,tempx1);
  xgrid((i-1)*norder+1:i*norder) = tempx;
  w((i-1)*norder+1:i*norder) = tempw;
end


%% evaluate function
n = numel(xgrid);
Ns = ones(1,d)*n;
ftt = amen_cross2(Ns, f, tol, xgrid);

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
Err = abs(S-Itrue);
fprintf('Err = %.50e \n',Err);
end


