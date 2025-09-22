function mom_for_node = cal_mom_for_node_3D(k,j,i,nnode_ext,ttmesh,vol_cell)
dim = 3;
mom_for_node = zeros(1,dim);
if (k >= 2 && k <= nnode_ext(3)-1) && ...
    (j >= 2 && j <= nnode_ext(2)-1) && ...
    (i >= 2 && i <= nnode_ext(1)-1)
  dmass = vol_cell * mean(full(ttmesh.rho_for_3dcell(k-1:k,j-1:j,i-1:i)),'all');
  mom_for_node(1:dim) = dmass * full(ttmesh.vel_for_3dnode(k, j, i, 1:dim));
end
end %function
