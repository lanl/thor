function R = calculateTimeDerivative(x,K,j1,j2,tint_order,dt)
  %
  R = zeros([size(x,1),size(x,2),size(x,3),size(x,4),K]);
  %
  j=j1:j2;
  %
  if(tint_order==1)
    R(:,:,:,:,1:K) = x(:,:,:,:,j) - x(:,:,:,:,j-1);
  elseif(tint_order==2)
    R(:,:,:,:,1:K) = (3/2)*x(:,:,:,:,j) - (2)*x(:,:,:,:,j-1) + (1/2)*x(:,:,:,:,j-2);
  elseif(tint_order==3)
    R(:,:,:,:,1:K) = (11/6)*x(:,:,:,:,j) - (3)*x(:,:,:,:,j-1) + (3/2)*x(:,:,:,:,j-2) - (1/3)*x(:,:,:,:,j-3);
  elseif(tint_order==4)
    R(:,:,:,:,1:K) = (25/12)*x(:,:,:,:,j) - (4)*x(:,:,:,:,j-1) + (3)*x(:,:,:,:,j-2) - (4/3)*x(:,:,:,:,j-3) + (1/4)*x(:,:,:,:,j-4);
  else
    error("tint_order must be less than or equal to 4");
  end
  %
  R(:) = R(:)/dt; 
  %
end