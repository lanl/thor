function error = check_tt_error(X,Xtt)
error = norm(X(:) - full(Xtt));
end