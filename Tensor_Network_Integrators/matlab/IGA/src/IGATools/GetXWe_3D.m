function [ctpxe, ctpye, ctpze, we] = GetXWe_3D(IndexE, ctpxv, ctpyv, ctpzv, wn)

ctpxe = ctpxv(IndexE);
ctpye = ctpyv(IndexE);
ctpze = ctpzv(IndexE);
we    = wn(IndexE);

end
