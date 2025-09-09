function S = pad_sigmatSF_1D_fn(S,imu,factor)
%1 - negative - pad the end (right/bottom)
%2 - positive - pad the beginning (left/top)

%%
temp = factor*ones(size(S(1,:,:,:)));
if imu==1
  % S = [S, factor*ones(nx,1)];
  S = cat(1,S,temp);
else
  S = cat(1,temp,S);
end

end %function