function [DMA,XPt] = compute_derivative_mat_fn(SN,type)

if strcmp(type,'Legendre')
  [XPt,VN] = zelegl(SN); % COMPUTES THE NODES RELATIVE
  XPt = XPt';
  DMA = dmlegl(SN,XPt,VN); % COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
elseif strcmp(type,'Chebyshev')
  XPt = ZECHGL(SN)';
  DMA = dmchgl(SN,XPt);
else
  error(['Wrong basis function type !!!\n' ...
    'Options are "Legendre" or "Chebyshev"'])
end