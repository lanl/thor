function Exc = tt_intp_Ec_fn(Ex, dim, tol)
% compute the cell center for field E in tensor train format
d1 = dim(1);
d2 = dim(2);

G = core2cell(Ex);
G{d1} = 0.5*(G{d1}(:,2:end,:) + G{d1}(:,1:end-1,:));
G{d2} = 0.5*(G{d2}(:,2:end,:) + G{d2}(:,1:end-1,:));
Exc = round(cell2core(tt_tensor,G),tol);


