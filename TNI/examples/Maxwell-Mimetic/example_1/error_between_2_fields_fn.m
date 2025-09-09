function [err,nrm] = error_between_2_fields_fn(X1, Y1, Z1, X2, Y2, Z2)
% compare the error between 2 fields

err = norm(X1(:)-X2(:)) + norm(Y1(:) - Y2(:)) + norm(Z1(:)-Z2(:));
nrm = norm(X1(:)) + norm(Y1(:)) + norm(Z1(:));

end