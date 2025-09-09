function S = pad_sigmat_1D_fn(S,imu,factor)
%1 - negative - pad the end (right/bottom)
%2 - positive - pad the beginning (left/top)

n = size(S,1);

if imu==1
  S = [S, factor*ones(1,n)];
else
  S = [factor*ones(1,n), S];
end

end %function