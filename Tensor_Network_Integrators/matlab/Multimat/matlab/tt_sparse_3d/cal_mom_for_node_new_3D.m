function mom_new = cal_mom_for_node_new_3D( ...
  k, j, i, vel_for_3dnode, mom_for_node, dir, dt, dx, small, dim, nnode_ext)
% CAL_MOM_FOR_NODE_NEW_3D
% Update momentum at node (k,j,i) by applying:
%   - its own outgoing flux, and
%   - incoming flux from BOTH neighboring nodes along 'dir' if they advect into (k,j,i).
%
% Inputs:
%   vel_for_3dnode(k,j,i,dir) : nodal velocity component along 'dir'
%   mom_for_node(k,j,i,idx)   : momentum components stored at nodes
%   dir ∈ {1,2,3}, dt, dx(1:3), small, dim, nnode_ext(1:3)
%
% Output:
%   mom_new(1,dim) : updated momentum vector at (k,j,i)

mom_new = zeros(1, dim);

% Helper to check node-index is an interior node (matches original loops: 2 .. N-1)
in_bounds = @(kk,jj,ii) ...
  kk >= 2 && kk <= nnode_ext(3)-1 && ...
  jj >= 2 && jj <= nnode_ext(2)-1 && ...
  ii >= 2 && ii <= nnode_ext(1)-1;

% Precompute outgoing fraction at (k,j,i)
dist_out = vel_for_3dnode(k, j, i, dir) * dt;
frac_out = abs(dist_out) / dx(dir);
has_out  = (frac_out >= small);

% Neighbor indices along 'dir'
km = k; jm = j; im = i;   % lower neighbor (−1 along dir)
kp = k; jp = j; ip = i;   % upper neighbor (+1 along dir)
switch dir
  case 1
    im = i - 1; ip = i + 1;
  case 2
    jm = j - 1; jp = j + 1;
  case 3
    km = k - 1; kp = k + 1;
end

% Loop momentum components
for idx = 1:dim
  delta = 0.0;

  % 1) Outgoing flux from (k,j,i)
  if has_out
    dmom_out = frac_out * mom_for_node(k, j, i, idx);
    delta = delta - dmom_out;
  end

  % 2) Incoming from LOWER neighbor (−1 along dir) if that neighbor advects toward +dir
  if in_bounds(km, jm, im)
    dist_m = vel_for_3dnode(km, jm, im, dir) * dt;
    if dist_m > 0.0     % lower neighbor sends to +dir -> into (k,j,i)
      frac_m = abs(dist_m) / dx(dir);
      if frac_m >= small
        dmom_m = frac_m * mom_for_node(km, jm, im, idx);
        delta = delta + dmom_m;
      end
    end
  end

  % 3) Incoming from UPPER neighbor (+1 along dir) if that neighbor advects toward −dir
  if in_bounds(kp, jp, ip)
    dist_p = vel_for_3dnode(kp, jp, ip, dir) * dt;
    if dist_p < 0.0     % upper neighbor sends to -dir -> into (k,j,i)
      frac_p = abs(dist_p) / dx(dir);
      if frac_p >= small
        dmom_p = frac_p * mom_for_node(kp, jp, ip, idx);
        delta = delta + dmom_p;
      end
    end
  end

  % Final updated momentum at (k,j,i)
  mom_new(idx) = mom_for_node(k, j, i, idx) + delta;
end
end


% function mom_new = cal_mom_for_node_new_3D(k, j, i, vel_for_3dnode, mom_for_node, dir,...
%   dt, dx, small, dim, nnode_ext)
% % CAL_MOM_FOR_NODE_NEW_3D
% % Computes updated momentum for node (k,j,i) by incorporating net momentum flux
% % from neighbors, but returns update only at that node.
%
%   mom_new = zeros(1, dim);  % Resulting momentum at (k,j,i)
%
%   for idx = 1:dim
%     delta_mom = 0.0;
%
%     % Outgoing flux from current node
%     dist_out = full(vel_for_3dnode(k, j, i, dir)) * dt;
%     frac_out = abs(dist_out) / dx(dir);
%
%     % Only apply outgoing momentum if valid
%     if frac_out >= small
%       dmom_out = frac_out * mom_for_node(k, j, i, idx);
%       delta_mom = delta_mom - dmom_out;
%
%       % Compute neighbor offset based on direction
%       kn = k; jn = j; in = i;
%       if dir == 1
%         in = i - sign(dist_out);
%       elseif dir == 2
%         jn = j - sign(dist_out);
%       elseif dir == 3
%         kn = k - sign(dist_out);
%       end
%
%       % Check bounds for neighbor before accessing
%       if kn >= 2 && kn <= nnode_ext(3)-1 && ...
%          jn >= 2 && jn <= nnode_ext(2)-1 && ...
%          in >= 2 && in <= nnode_ext(1)-1
%
%         dist_in = full(vel_for_3dnode(kn, jn, in, dir)) * dt;
%         frac_in = abs(dist_in) / dx(dir);
%
%         if frac_in >= small
%           dmom_in = frac_in * mom_for_node(kn, jn, in, idx);
%           delta_mom = delta_mom + dmom_in;
%         end
%       end
%     end
%
%     % Final momentum update
%     mom_new(idx) = mom_for_node(k, j, i, idx) + delta_mom;
%   end
% end
% % function mom_new = cal_mom_for_node_new_3D(k, j, i, vel_for_3dnode, mom_for_node, dir,...
% %   dt, dx, small, dim, nnode_ext)
% % % CAL_MOM_FOR_NODE_NEW_3D
% % % Computes updated momentum for node (k,j,i) by incorporating net momentum flux
% % % from neighbors, but returns update only at that node.
% %
% %   mom_new = zeros(1, dim);  % Resulting momentum at (k,j,i)
% %
% %   for idx = 1:dim
% %     delta_mom = 0.0;
% %
% %     % Outgoing flux from current node
% %     dist_out = full(vel_for_3dnode(k, j, i, dir)) * dt;
% %     frac_out = abs(dist_out) / dx(dir);
% %     if frac_out >= small
% %       dmom_out = frac_out * mom_for_node(k, j, i, idx);
% %       delta_mom = delta_mom - dmom_out;
% %     end
% %
% %     % Incoming flux from neighbor in opposite direction
% %     kn = k; jn = j; in = i;
% %     if dir == 1, in = i - sign(dist_out); end
% %     if dir == 2, jn = j - sign(dist_out); end
% %     if dir == 3, kn = k - sign(dist_out); end
% %
% %     % Check neighbor is in bounds
% %     if kn >= 3 && kn <= nnode_ext(3)-2 && ...
% %        jn >= 3 && jn <= nnode_ext(2)-2 && ...
% %        in >= 3 && in <= nnode_ext(1)-2
% %       dist_in = full(vel_for_3dnode(kn, jn, in, dir)) * dt;
% %       frac_in = abs(dist_in) / dx(dir);
% %       if frac_in >= small
% %         dmom_in = frac_in * mom_for_node(kn, jn, in, idx);
% %         delta_mom = delta_mom + dmom_in;
% %       end
% %     end
% %
% %     mom_new(idx) = mom_for_node(k, j, i, idx) + delta_mom;
% %   end
% % end
