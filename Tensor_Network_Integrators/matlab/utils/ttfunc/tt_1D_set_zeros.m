function X = tt_1D_set_zeros(X, dim, idx)
% this function will set corresponding slices in idx in the dimension/core 'dim'
% to zeros

G = X{dim};
G(:,idx,:) = 0;
X{dim} = G;