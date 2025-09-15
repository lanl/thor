function X = expand3d(X,d)
%d = [dim ix nx]
dim = d(1);
ix= d(2);
nx = d(3);

%get the core in the dim dimension
r = X.r;
Gi = zeros(r(dim),nx,r(dim+1));
Gi(:,ix,:) = X{dim};
X{dim} = Gi;% reassign the expanded core
% %% old code 

% G = core2cell(X);
% %for now, this is only for the first and second dimension
% v = zeros(r(dim),nx,r(dim+1));
% v(:,ix,:) = G{dim};
% if dim==1
%   X = cell2core(tt_tensor,{v,G{2},G{3}});
% elseif dim==2
%   X = cell2core(tt_tensor,{G{1},v,G{3}});
% elseif dim==3
%   X = cell2core(tt_tensor, {G{1},G{2},v});
% end

end