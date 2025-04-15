function [Qstencil,idx_m] = getWeno5StencilNeq(Q,d)
  %
  Neq = length(Q);
  %
  Qstencil = cell(1,5*Neq);
  %
  rmax  = 0;
  idx_m = 0;
  %
  for eq=1:Neq
    %
    for j=1:5
      %
      idx = eq + Neq*(j-1);
      %
      Qstencil{idx} = applyShift(Q{eq},d,j-3);
      %
      if(Qstencil{idx}.r>rmax)
        %
        rmax  = Qstencil{idx}.r;
        %
        idx_m = idx;
        %
      end
      %
    end
    %
  end
  %