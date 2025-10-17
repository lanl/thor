function Y = FD_tt_ops(X,ops,Cidx,ops_factor)

%apply Finite Difference operators on X in TT format
%ops is a cell array of operations -string format. Avaliable ops are:
%'diff','intp','none'
%Cidx is a cell array containing sets of indices for each operation. In
%each cell are two set of idxs for that corresponding operation in ops.

if isa(X,'tt_tensor')
  dim = X.d;
  G = core2cell(X);
elseif isa(X,'cell')
  dim = numel(X);
  G = X;
end
for i = 1:dim
  I1 = Cidx{i}{1};
  I2 = Cidx{i}{2};
  switch ops{i}
    case 'diff'
      G{i} = G{i}(:,I2,:) - G{i}(:,I1,:);
    case 'intp'
      G{i} = G{i}(:,I2,:) + G{i}(:,I1,:);
    case 'none'
      G{i} = G{i}(:,I1,:);
    otherwise
      error("Unexpected FD operation. Are you kidding me !!!")
  end
  G{i} = ops_factor(i)*G{i};
end

Y = cell2core(tt_tensor,G);