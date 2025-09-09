function err = tt_sampling_error(X,Ytt,npts)
% randomly sample npts points in X and Y to check for the error
%X can be a full matrix or TT, Y is a TT
err = 0;

if isa(Ytt,'tt_matrix')
  N = prod(size(Ytt),1);
  n1 = N(1); n2 = N(2);

else
  n1 = prod(size(Ytt));
  n2 = 1;
end

if issparse(X)
  [J1,J2] = find(X);
  for i = 1:min(npts,numel(J1))
    j1 = J1(i);
    j2 = J2(i);
    err = err + abs(X(j1,j2) - Ytt(j1,j2));
  end
  err = err/min(npts,numel(J1));
else
  I = randi(n1*n2,npts,1); %sample indices

  for i = 1:npts
    [j1,j2] = ind2sub([n1,n2],I(i));
    err = err + abs(X(j1,j2) - Ytt(j1,j2));
  end
  err = err/npts;
end


end