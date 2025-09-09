function tt=create_tt(n,r,d)
  %
  tt     = tt_tensor;
  %
  tt.n      = n;
  tt.d      = d;
  tt.r      = r;
  tt.ps     = cumsum([1;n.*r(1:d).*r(2:d+1)]);
  tt.core   = zeros(tt.ps(d+1)-1,1);
  %
end