function [CC, k] = KnotInsertMC(knot,t)
 [m, ~] = EvaluateKnot(knot);
 n      = length(t);
 tt     = sort(t);
 k      = knot;
 CC     = eye(m);
 
 for i=1:n
  CCm = KnotInsertC(k, tt(i));
  CC = CCm*CC;
  k(length(k)+1) = tt(i);
  k = sort(k);
 end
 
end

