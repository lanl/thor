function [G,r] = tt_orthogonolize(G,n,r,type)
  %
  % type =  1: Left to Right
  % type = -1: Right to Left
  %
  d = numel(G);
  %
  if(type==1)
    %
    for i=1:d
      %
      cr = reshape(G{i},[r(i)*n(i),r(i+1)]);
      %
      [tempq,tempr] = qr(cr,"econ");
      %
      G{i} = reshape(tempq,r(i),n(i),r(i+1));
      %
      if(i<d)
        G{i+1} = tensorprod(tempr,G{i+1},2,1);
      else
        r = tempr;
      end
      %
    end
    %
  elseif(type==-1)
    %
    for i=d:-1:1
      %
      cr = reshape(G{i},[r(i),n(i)*r(i+1)]);
      %
      [tempq,tempr] = qr(cr',"econ");
      %
      tempq = tempq';
      tempr = tempr';
      %
      G{i} = reshape(tempq,r(i),n(i),r(i+1));
      %
      if(i>1)
        G{i-1} = tensorprod(G{i-1},tempr,3,1);
      else
        r = tempr;
      end
      %
    end
    %
  else
    %
    error("Invalid value: type = %d",type);
    %
  end
  %
end