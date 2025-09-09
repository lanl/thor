function tt_assert_error(X,Xtt,tol,txt)
if nargin<4
  txt = " checking error ";
end
if isa(Xtt,'tt_matrix')
  error = norm(X - full(Xtt));
else
  error = norm(X(:) - full(Xtt));
end
if tol>0
  assert(error<tol,sprintf('Error in %s= %.5e \n', txt, error))
end
printTextWithOK(sprintf('%s-err=%.2e',txt,error),55);