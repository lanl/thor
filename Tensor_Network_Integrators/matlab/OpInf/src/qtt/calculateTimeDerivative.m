function R = calculateTimeDerivative(x,K,j1,j2,tint_order,dt,eps_tt,isQtt)
  %
  R = x;
  %
  G      = core2cell(R);
  G{end} = G{end}(:,1:K,:);
  %
  cr = x{x.d};
  %
  j=j1:j2;
  %
  if(tint_order==1)
    G{end} = cr(:,j,:) - cr(:,j-1,:);
  elseif(tint_order==2)
    G{end} = (3/2)*cr(:,j,:) - (2)*cr(:,j-1,:) + (1/2)*cr(:,j-2,:);
  elseif(tint_order==3)
    G{end} = (11/6)*cr(:,j,:) - (3)*cr(:,j-1,:) + (3/2)*cr(:,j-2,:) - (1/3)*cr(:,j-3,:);
  elseif(tint_order==4)
    G{end} = (25/12)*cr(:,j,:) - (4)*cr(:,j-1,:) + (3)*cr(:,j-2,:) - (4/3)*cr(:,j-3,:) + (1/4)*cr(:,j-4,:);
  else
    error("tint_order must be less than or equal to 4");
  end
  %
  G{end} =G{end}/dt; 
  %
  R = cell2core(tt_tensor,G);
  %
  if(isQtt)
    %
    n = [R.n(1:end-2)' R.n(end-1) factor(R.n(end))];
    %
    R = tt_reshape(R,n);
    %
  end
  %
  R = round(R,eps_tt);
  %
end