function error = check_tt_rel_error(X,Xtt)
error = norm(X(:) - full(Xtt))/norm(X(:));
end