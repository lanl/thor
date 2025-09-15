function S = pad_sigmat_fn(S,imu,jeta,factor)
%1 - negative - pad the end (right/bottom)
%2 - positive - pad the beginning (left/top)

n = size(S,1);
if jeta==1
  S = [S, factor*ones(n,1)];
else
  S = [factor*ones(n,1) S];
end

if imu==1
  S = [S; factor*ones(1,n+1)];
else
  S = [factor*ones(1,n+1);S];
end

end %function