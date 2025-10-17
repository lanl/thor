function X = take_out_small_vof_3D(k, j, i, obj, mat, nbdry, ncell_bdry, nmat, vol_cell)

vol_for_cell = full(obj.vol_for_cell(k,j,i,:))';
mass_for_cell = full(obj.mass_for_cell(k,j,i,:))';
ener_for_cell = full(obj.ener_for_cell(k,j,i,:))';

if (k >= nbdry+1 && k <= ncell_bdry(3)+1) && ...
    (j >= nbdry+1 && j <= ncell_bdry(2)+1) && ...
    (i >= nbdry+1 && i <= ncell_bdry(1)+1)
  vol_sum = 0.0;

  % Remove small volume fractions below threshold
  for m = 1:nmat
    frac = vol_for_cell(m) / vol_cell;
    if frac < mat.vfmin
      vol_for_cell(m) = 0.0;
      mass_for_cell(m) = 0.0;
      ener_for_cell(m) = 0.0;
    else
      vol_sum = vol_sum + obj.vol_for_cell(k, j, i, m);
    end
  end

  % Normalize remaining materials to maintain conservation
  if vol_sum > 0  % Avoid division by zero
    frac = vol_cell / vol_sum;
    vol_for_cell(1:nmat) = ...
      vol_for_cell(1:nmat) * frac;
  end
end %if

X = [vol_for_cell, mass_for_cell, ener_for_cell];


