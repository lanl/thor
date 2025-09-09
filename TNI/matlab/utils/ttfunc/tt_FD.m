function X = tt_FD(X,dim)
% 12/14/21 implement the finite difference on tensor train X in the dimension
% 'dim'. ONLY SHIFT THE CORRESPONDING CORE

%get the core in that dimension
Gi = X{dim};
% compute finite difference on the core
Gi = Gi(:,2:end,:) - Gi(:,1:end-1,:);
X{dim} = Gi;
    
end