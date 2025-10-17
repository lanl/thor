function R = calculateTimeDerivative(x,K,n,p,s,j1,j2,isA,isF,tint_order,dt,g1,g2)
  %
  Rsz = 0;
  %
  if(g1>0)
    if(isA)
      Rsz = Rsz + n;
    end
    Rsz = Rsz + p;
  end
  %
  if(g2>0 && isF) 
    Rsz = Rsz + s;
  end
  %
  R = zeros(K+Rsz,size(x,2));
  %
  j=j1:j2;
  %
  if(tint_order==1)
    R(1:K,:) = x(j,:) - x(j-1,:);
  elseif(tint_order==2)
    R(1:K,:) = (3/2)*x(j,:) - (2)*x(j-1,:) + (1/2)*x(j-2,:);
  elseif(tint_order==3)
    R(1:K,:) = (11/6)*x(j,:) - (3)*x(j-1,:) + (3/2)*x(j-2,:) - (1/3)*x(j-3,:);
  elseif(tint_order==4)
    R(1:K,:) = (25/12)*x(j,:) - (4)*x(j-1,:) + (3)*x(j-2,:) - (4/3)*x(j-3,:) + (1/4)*x(j-4,:);
  else
    error("tint_order must be less than or equal to 4");
  end
  %
  R(:) = R(:)/dt; 
  %
end