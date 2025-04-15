function X = tt_zeropad_fn(X,dims,nrows)
% this function will zero pad X on the dim in 'dims'
% numrows: dimension dimsx2. 
% in each row, [numrow padded in the beginning, numrow padded in the end]

for i = 1:numel(dims)
  d = dims(i);
  r = nrows(i,:);
  if d == ndims(X)
    if r(1)>0
      X{d} = cat(2,zeros(size(X{d},1),r(1)), X{d});
    end
    if r(2)>0
      X{d} = cat(2,X{d},zeros(size(X{d},1),r(2)));
    end
  else
    if r(1) > 0
      X{d} = cat(2,zeros(size(X{d},1),r(1),size(X{d},3)), X{d});
    end
    if r(2) > 0
      X{d} = cat(2, X{d}, zeros(size(X{d},1),r(2),size(X{d},3)));
    end
    
  end
    
end