function RHS=calculateReducedRHS(t,x,u,A,F,B,isA,isF)
  %
  n = length(x);
  %
  if(~isempty(B))
    RHS = B*u(t);
  else
    RHS = zeros(n,1);
  end 
  %
  if(isA)
    RHS(:) = RHS + A*x; 
  end
  %
  if(isF)
    %
    x2 = getx2(x,n);
    %
    RHS(:) = RHS + F*x2;
    %
  end
  %
end
%
function x2 = getx2(x,n)
  %
  x2 = zeros(n*(n+1)/2,1);
  %
  i1 = 1;
  %
  for i=1:n
    %
    i2 = i1 + i - 1;
    %
    x2(i1:i2) = x(i)*x(1:i);
    %
    i1 = i2 + 1;
    %
  end
  %
end
%