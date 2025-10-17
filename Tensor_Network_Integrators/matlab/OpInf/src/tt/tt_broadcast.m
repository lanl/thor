function tt = tt_broadcast(tt, dims, N)
  %
  if(length(dims)~=length(N))
    error("length(dims) must be the same as length(N)");
  end
  %
  G = core2cell(tt);
  %
  for idx=1:length(dims)
    %
    i = dims(idx);
    n = N(idx);
    %
    if(size(G{i},2)~=1)
      error("broadcasting is allowed only when n=1");
    end
    %
    G{i} = G{i}(:,ones(1,n),:);
    %
  end
  %
  tt = cell2core(tt_tensor,G);
  %
end