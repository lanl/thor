function Xpadded = tt_pad_zeros(X)

dim = X.d;
n = X.n;
r = X.r;
G = core2cell(X); % get cores from X
Gnew = cell(dim,1); % list of new cores

for i = 1:dim
  Gnew{i} = zeros(r(i),n(i)+2,r(i+1));
  Gnew{i}(:,2:n(i)+1,:) = G{i};
end

Xpadded = cell2core(tt_tensor, Gnew);
