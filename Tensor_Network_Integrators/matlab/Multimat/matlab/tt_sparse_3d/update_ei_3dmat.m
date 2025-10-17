function ei_3dmatnew = update_ei_3dmat(k,j,i,m,ei_3dmat,obj,...
  nbdry,ncell_bdry,tiny)

% Check bounds
if (k < nbdry + 1 || k > ncell_bdry(3) + 1 || ...
    j < nbdry + 1 || j > ncell_bdry(2) + 1 || ...
    i < nbdry + 1 || i > ncell_bdry(1) + 1)
  ei_3dmatnew = ei_3dmat(k,j,i,m);
  return;  % Skip if (k,j,i) is out of valid interior region

elseif obj.ener_for_cell(k, j, i, m)<1e-13
    ei_3dmatnew = 0;
  return; 
else
  ei_3dmatnew = obj.ener_for_cell(k, j, i, m) / ...
    (tiny + obj.vol_for_cell(k, j, i, m));
end
end