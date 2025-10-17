function vf_3dmatnew = update_vf_3dmat(k,j,i,m,vf_3dmat,obj,vol_cell,...
  nbdry,ncell_bdry)

% Check bounds
if (k < nbdry + 1 || k > ncell_bdry(3) + 1 || ...
    j < nbdry + 1 || j > ncell_bdry(2) + 1 || ...
    i < nbdry + 1 || i > ncell_bdry(1) + 1)
  vf_3dmatnew = vf_3dmat(k,j,i,m);
  return;  % Skip if (k,j,i) is out of valid interior region
else
  vf_3dmatnew = min(1.0, obj.vol_for_cell(k, j, i, m) / vol_cell);
end
end