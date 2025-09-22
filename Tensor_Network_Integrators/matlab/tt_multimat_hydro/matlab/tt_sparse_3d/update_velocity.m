function vel_for_3dnode = update_velocity(k,j,i,nbdry,nnode_bdry,ttmesh,...
  mom_for_3dnode_newtt,rho_average,vol_cell,dim,small)

if (k >= nbdry && k <= nnode_bdry(3)) && ...
    (j >= nbdry && j <= nnode_bdry(2)) && ...
    (i >= nbdry && i <= nnode_bdry(1))
  rho_sum = mean(full(ttmesh.rho_for_3dcell(k-1:k, j-1:j, i-1:i)), 'all');
  if rho_sum / rho_average < small
    vel_for_3dnode(1:dim) = 0.0;
  else
    % Update velocity based on momentum and density
    vel_for_3dnode(1:dim) = ...
      full(mom_for_3dnode_newtt(k, j, i, 1:dim))./(vol_cell * rho_sum);
  end
else 
  vel_for_3dnode(1:dim) = full(ttmesh.vel_for_3dnode(k,j,i,1:dim));
end
end
