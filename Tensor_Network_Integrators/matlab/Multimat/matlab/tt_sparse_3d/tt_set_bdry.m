function Xnew = tt_set_bdry(X,dim,Idx)

% Set layer i1 using layer i0's value, for dimension dim

i0 = Idx(1);
i1 = Idx(2);
G = core2cell(X);
G{dim}(:,i1,:) = G{dim}(:,i0,:);

Xnew = cell2core(tt_tensor,G);

end