function error = check_tt_error(X,Xtt,txt)

if isa(Xtt,'tt_matrix')
  if isa(X,'tt_matrix')
    error = norm(full(X) - full(Xtt));
  else
    error = norm(X - full(Xtt));
  end
else
  if isa(X,'tt_tensor')
    error = norm(full(X) - full(Xtt));
  else
    error = norm(X(:) - full(Xtt));
  end
end

if nargin>2
  fprintf('Error in %s = %.5e \n',txt,error);
end