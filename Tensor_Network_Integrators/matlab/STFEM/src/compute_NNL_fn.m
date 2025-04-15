function NNL = compute_NNL_fn(Ex)
Nx = Ex + 1;
NNL= sparse(Nx, 2*Ex);
for ii=2:Ex
  NNL(ii,ii+(ii-2))=1; NNL(ii,ii+(ii-2)+1)=1;
end
NNL(1,1)=1; NNL(Nx, 2*Ex)=1;
