function angttPsi1 = compute_angflux_fn(ttPsi1)

G = core2cell(ttPsi1);
G{5} = sum(G{5},2);
G{4} = tensorprod(G{4},G{5},3,1);
angttPsi1 = cell2core(tt_tensor,G(1:4));

end