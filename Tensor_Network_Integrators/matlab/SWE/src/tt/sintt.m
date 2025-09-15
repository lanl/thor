function tt=sintt(x,i)
  %
  tt = tt_ones(x.n,x.d);
  %
  i1 = tt.ps(i);
  i2 = tt.ps(i+1)-1;
  %
  tt.core(i1:i2) = sin(x.core(i1:i2));
  %
end