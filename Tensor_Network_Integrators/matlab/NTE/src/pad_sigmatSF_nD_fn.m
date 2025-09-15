function S = pad_sigmatSF_nD_fn(S,imu,jeta,factor)
%1 - negative - pad the end (right/bottom)
%2 - positive - pad the beginning (left/top)

temp = factor*ones(size(S(:,1,:,:)));
if jeta==1
  % S = [S, factor*ones(nx,1)];
  S = cat(2,S,temp);
else
  S = cat(2,temp,S);
end

%%
temp = factor*ones(size(S(1,:,:,:)));
if imu==1
  % S = [S, factor*ones(nx,1)];
  S = cat(1,S,temp);
else
  S = cat(1,temp,S);
end

end %function