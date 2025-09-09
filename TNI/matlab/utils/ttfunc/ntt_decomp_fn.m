function [Xtt] = ntt_decomp_fn(X,r)
% nonnegative tensor train decomposition with rank r
% return a TT-tensor object

N = size(X);
d = numel(N);
Ks = r;
node = Node(char(64 + [1:d]));
kstr = char(64 + numel(N) + [1:d-1]);

counter = numel(Ks);
node = tensor_train_tree_fn(node, Ks, kstr, counter);

%% pass the tree into recursive_tensor_decomp
node = recursive_tensor_decomp(@decompfunc,double(X),node);

%% construct the TT tensor from node
G = cell(d,1);
crnode = node;
for i = 1:d-1
  G{i} = crnode.left.factor;
  crnode = crnode.right;
end
G{d} = crnode.factor;
% reshape G{1}
sz = [1, size(G{1})];
G{1} = reshape(G{1},sz);
Xtt = cell2core(tt_tensor, G);

%get err

end

%%
function [W,H] = decompfunc(X,W,H)
[W,H] = nmf_MU_Fro(X,W,H,0,5000);
end
