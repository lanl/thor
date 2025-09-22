function error = check_tt_rel_errorh(X, Xtt, h)
Xttf = full(Xtt);
    diff = X(:) - Xttf(:);
    error = sqrt(h * sum(diff.^2)) / sqrt(h * sum(X(:).^2));
end
