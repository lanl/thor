function pres_3dmatnew = update_pres_3dmat(k,j,i,m,mesh,...
  nbdry,ncell_bdry,is_solid,gamma_ea_mat)

% Check bounds
if (k < nbdry + 1 || k > ncell_bdry(3) + 1 || ...
    j < nbdry + 1 || j > ncell_bdry(2) + 1 || ...
    i < nbdry + 1 || i > ncell_bdry(1) + 1)
  pres_3dmatnew = mesh.pres_3dmat(k,j,i,m);
  return;  % Skip if (k,j,i) is out of valid interior region
else
  if is_solid(m)
    pres_3dmatnew = ...
      eos.p_mie_gruneisen(mesh.rho_3dmat(k, j, i, m), ...
      mesh.ei_3dmat(k, j, i, m));
  else
    pres_3dmatnew = ...
      (gamma_ea_mat(m) - 1.0) * mesh.ei_3dmat(k, j, i, m);
  end
end
end