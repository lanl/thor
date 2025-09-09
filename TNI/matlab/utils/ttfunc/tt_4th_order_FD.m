function X = tt_4th_order_FD(X,dim)
% 12/14/21 implement the finite difference on tensor train X in the dimension
% 'dim'. ONLY SHIFT THE CORRESPONDING CORE

%get the core in that dimension
Gi = X{dim};

% compute finite difference on the core
% Gi = Gi(:,2:end,:) - Gi(:,1:end-1,:);
temp = zeros(X.r(dim),X.n(dim)-1, X.r(dim+1));
temp(:,2:end-1,:) = 1/48*Gi(:,1:end-3,:) - 9/16*Gi(:,2:end-2,:)... 
+ 9/16*Gi(:,3:end-1,:) - 1/48*Gi(:,4:end,:);
% pad zeros
% temp = [zeros]
X{dim} = temp;
    
end