function [V,n,retained_energy] = calculatePOD(x,param)
  %
  [U,S,~] = svd(x,"econ");
  %
  S = diag(S);
  %
  POD_energy = cumsum(S.^2) / sum(S.^2);
  %
  if(param<1)
    %
    retained_energy = param;
    %
    n = find(POD_energy > retained_energy, 1);
    %
  else
    n = min(length(POD_energy),param);
  end
  %
  retained_energy = POD_energy(n);
  %
  V = U(:,1:n);
  %
end