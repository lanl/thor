function Att = tt_append_operator_fn(Att,Bmat)
% Att is in the TT-matrix format
% Bmat is a matrix that will be appended in the end of Att
G = core2cell(Att);
G = [G; {reshape(Bmat,[1,size(Bmat)])}];

Att = cell2core(tt_matrix,G);
end