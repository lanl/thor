function [ctpxv, ctpyv, ctpzv, IndexE, ctpxe, ctpye, ctpze, we, m3D] = ...
                                       IGAinfo(ctpxn, ctpyn, ctpzn, wn,...
                                               knot1n, knot2n, knot3n)
% Vectorize control points
[ctpxv, ctpyv, ctpzv] = Vectorizectp3D(ctpxn, ctpyn, ctpzn);

% Conectivity and support domain of each element
%[IndexE, m3D] = N3DIndex(knot1n, knot2n, knot3n);
[IndexE, m3D] = N3DIndex_optimized(knot1n, knot2n, knot3n);

% Control points and weigth factors of each element
[ctpxe, ctpye, ctpze, we] = GetXWe_3D(IndexE, ctpxv, ctpyv, ctpzv, wn);
end