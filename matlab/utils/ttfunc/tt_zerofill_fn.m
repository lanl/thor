function X = tt_zerofill_fn(X,idxcell)
% idxcell is a cell array containing indices for non-zeros areas. It could
% be a number or an array
% This function will set other elements to zeros.

%check if the input is valid
Nd = ndims(X);
if ndims(X) ~= numel(idxcell)
  error('the number of index cells are incorrect');
end

G = core2cell(X);

for i = 1: Nd
  In = idxcell{i};
  I0 = setdiff([1:X.n(i)],In);
  if i == Nd
    G{i}(:,I0) = 0;
  else
    G{i}(:,I0,:) = 0;
  end
end

X = cell2core(X,G);

