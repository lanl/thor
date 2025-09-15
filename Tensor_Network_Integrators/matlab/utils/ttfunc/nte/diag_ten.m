function T = diag_ten(r,n)
%make a diagonal tensor rxnxr where each slice is an identity
T = zeros(r,n,r);
for i = 1:n
  T(:,i,:) = eye(r);
end