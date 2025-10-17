function [x,t]=predict(x,u,A,F,B,m,Tfinal,dt,tint_order,isA,isF)
  %
  x = x(1:m,:);
  %
  if(isA)
    A = A(1:m,1:m);
  end
  %
  if(~isempty(B))
    B = B(1:m,:);
  end
  %
  if(isF)
    F = F(1:m,1:m*(m+1)/2);
  end
  %
  x = x(:,end);
  %
  if(tint_order<3)
    [t, x] = ode23(@(tt,xx)calculateReducedRHS(tt,xx,u,A,F,B,isA,isF), 0:dt:Tfinal, x);
    x=x';
  else
    [t, x] = ode45(@(tt,xx)calculateReducedRHS(tt,xx,u,A,F,B,isA,isF), 0:dt:Tfinal, x);
    x=x';
  end
  %
end