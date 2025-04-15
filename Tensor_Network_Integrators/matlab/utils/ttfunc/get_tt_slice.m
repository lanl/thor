function D = get_tt_slice(X,d)

%d =[dim ix iloc]
% in the 'dim' dimension, get ix slice and place it to iloc location

n = d(1);
ix = d(2);
iloc = d(3);

G = core2cell(X);

v = zeros(X.r(n),X.n(n),X.r(n+1));
v(:,iloc,:) = X{n}(:,ix,:);
G{n} = v;
D = cell2core(tt_tensor,G);
end