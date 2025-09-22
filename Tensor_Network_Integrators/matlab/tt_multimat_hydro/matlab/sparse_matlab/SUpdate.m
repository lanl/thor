classdef SUpdate < handle
  %C_UPDATE Class for updating mesh properties
  %   This class contains methods to compute and update various physical
  %   quantities on a mesh, such as divergence of velocity, viscous forces,
  %   and energy, among others.

  properties
    % Force arrays
    force_for_2dnode  % 2D force at nodes, initially an empty 2D array
    force_for_3dnode  % 3D force at nodes, initially an empty 3D array

    % Viscous force (qvis) arrays
    qvis_for_2dcell   % 2D viscous force at cells, initially an empty 2D array
    qvis_for_3dcell   % 3D viscous force at cells, initially an empty 3D array

    % Divergence of velocity (divu) arrays
    divu_for_2dcell   % 2D divergence of velocity at cells, initially an empty 2D array
    divu_for_3dcell   % 3D divergence of velocity at cells, initially an empty 3D array
  end

  methods

    function compute_divu(obj, mesh)
      % Compute the divergence of the velocity field
      % Inputs:
      %   - obj: instance of c_update
      %   - mesh: instance of the mesh class with properties defining the problem

      % Shortcuts
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      dx = mesh.dx(1:dim);
      ncell_ext = ncell + 2 * nbdry;  % Extended cell count including boundaries
      volinv = 1.0 / prod(dx(1:dim));  % Inverse of the cell volume

      if dim == 2
        % Initialize the divu_for_2dcell array if it is empty
        obj.divu_for_2dcell = zeros(ncell_ext(2), ncell_ext(1));
        % Calculate divergence for 2D
        for j = 1:ncell_ext(2)
          j1 = j + 1;
          for i = 1:ncell_ext(1)
            i1 = i + 1;

            dvx = 0.5 * ((mesh.vel_for_2dnode(j, i1, 1) + mesh.vel_for_2dnode(j1, i1, 1)) - ...
              (mesh.vel_for_2dnode(j, i, 1) + mesh.vel_for_2dnode(j1, i, 1)));

            dvy = 0.5 * ((mesh.vel_for_2dnode(j1, i, 2) + mesh.vel_for_2dnode(j1, i1, 2)) - ...
              (mesh.vel_for_2dnode(j, i, 2) + mesh.vel_for_2dnode(j, i1, 2)));

            obj.divu_for_2dcell(j, i) = volinv * (dx(2) * dvx + dx(1) * dvy);
          end
        end

      elseif dim == 3
        % Initialize the divu_for_3dcell array if it is empty
        obj.divu_for_3dcell = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));

        % Calculate divergence for 3D
        for k = 1:ncell_ext(3)
          k1 = k + 1;
          for j = 1:ncell_ext(2)
            j1 = j + 1;
            for i = 1:ncell_ext(1)
              i1 = i + 1;

              dvx = 0.25 * ((mesh.vel_for_3dnode(k, j, i1, 1) + mesh.vel_for_3dnode(k, j1, i1, 1) + ...
                mesh.vel_for_3dnode(k1, j, i1, 1) + mesh.vel_for_3dnode(k1, j1, i1, 1)) - ...
                (mesh.vel_for_3dnode(k, j, i, 1) + mesh.vel_for_3dnode(k, j1, i, 1) + ...
                mesh.vel_for_3dnode(k1, j, i, 1) + mesh.vel_for_3dnode(k1, j1, i, 1)));

              % dvy = 0.25 * ((mesh.vel_for_3dnode(k, j1, i, 2) + mesh.vel_for_3dnode(k, j1, i1, 2) + ...
              %   mesh.vel_for_3dnode(k1, j1, i, 2) + mesh.vel_for_3dnode(k1, j1, i1, 2)) - ...
              %   (mesh.vel_for_3dnode(k, j, i, 2) + mesh.vel_for_3dnode(k, j, i1, 2) + ...
              %   mesh.vel_for_3dnode(k1, j, i, 2) + mesh.vel_for_3dnode(k1, j1, i1, 2)));

              dvy = 0.25 * ( ...
                (mesh.vel_for_3dnode(k,  j1, i,  2) + mesh.vel_for_3dnode(k,  j1, i1, 2) + ...
                mesh.vel_for_3dnode(k1, j1, i,  2) + mesh.vel_for_3dnode(k1, j1, i1, 2)) - ...
                (mesh.vel_for_3dnode(k,  j,  i,  2) + mesh.vel_for_3dnode(k,  j,  i1, 2) + ...
                mesh.vel_for_3dnode(k1, j,  i,  2) + mesh.vel_for_3dnode(k1, j,  i1, 2)) );
              
              dvz = 0.25 * ((mesh.vel_for_3dnode(k1, j, i, 3) + mesh.vel_for_3dnode(k1, j, i1, 3) + ...
                mesh.vel_for_3dnode(k1, j1, i, 3) + mesh.vel_for_3dnode(k1, j1, i1, 3)) - ...
                (mesh.vel_for_3dnode(k, j, i, 3) + mesh.vel_for_3dnode(k, j, i1, 3) + ...
                mesh.vel_for_3dnode(k, j1, i, 3) + mesh.vel_for_3dnode(k, j1, i1, 3)));

              obj.divu_for_3dcell(k, j, i) = volinv*(dx(2)*dx(3)*dvx + dx(1)*dx(3)*dvy + dx(1)*dx(2)* dvz);
            end
          end
        end
      end
    end

    function compute_qvis(obj, mesh)
      % Compute the artificial viscosity (qvis) for the mesh
      % Inputs:
      %   - obj: instance of c_update
      %   - mesh: instance of the mesh class with properties defining the problem

      % Shortcuts
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      dx = mesh.dx(1);
      ncell_ext = ncell + 2 * nbdry;

      % Viscosity parameters
      c1 = 1.0;
      c2 = 0.1;

      if dim == 2
        % Initialize the qvis_for_2dcell array if it is empty
        obj.qvis_for_2dcell = zeros(ncell_ext(2), ncell_ext(1));

        % Calculate qvis for 2D cells
        for j = 1:ncell_ext(2)
          for i = 1:ncell_ext(1)
            if obj.divu_for_2dcell(j, i) < 0.0
              dv = dx * obj.divu_for_2dcell(j, i);
              obj.qvis_for_2dcell(j, i) = abs(dv) * mesh.rho_for_2dcell(j, i) * ...
                (c1 * mesh.cs_for_2dcell(j, i) - c2 * dv);
            end
          end
        end

      elseif dim == 3
        % Initialize the qvis_for_3dcell array if it is empty
        obj.qvis_for_3dcell = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));
        
        % Calculate qvis for 3D cells
        for k = 1:ncell_ext(3)
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              if obj.divu_for_3dcell(k, j, i) < 0.0
                dv = dx * obj.divu_for_3dcell(k, j, i);
                obj.qvis_for_3dcell(k, j, i) = abs(dv) * mesh.rho_for_3dcell(k, j, i) * ...
                  (c1 * mesh.cs_for_3dcell(k, j, i) - c2 * dv);
              end
            end
          end
        end
      end
    end

    %%
    % function compute_force(obj, mesh, direction)
    %   % Compute the force for the mesh
    %   % Inputs:
    %   %   - obj: instance of c_update
    %   %   - mesh: instance of the mesh class with properties defining the problem
    %   %   - direction: direction of force calculation (0, 1, or 2)
    %
    %   global tiny % global variable defined in main
    %
    %   % Extract necessary properties from the mesh
    %   dim = mesh.dim_prob;
    %   ncell = mesh.ncell_prob;
    %   nbdry = mesh.nbdry_prob;
    %   dx = mesh.dx(1:dim);
    %
    %   nnode_ext = ncell + 2 * nbdry + 1;  % Extended number of nodes
    %   vol_of_node = prod(0.5 * dx(1:dim));  % Volume associated with each node
    %   area_of_corner = prod(0.5 * dx(setdiff(1:dim, direction)));  % Area of corner perpendicular to the direction
    %
    %   % Check dimensionality and initialize force arrays
    %   if dim == 2
    %     % Initialize the force_for_2dnode array if it is empty
    %     if isempty(obj.force_for_2dnode)
    %       obj.force_for_2dnode = zeros(nnode_ext(2), nnode_ext(1));
    %     end
    %
    %     % Calculate force for 2D nodes
    %     for j = 2:nnode_ext(2) - 1
    %       j0 = j - 1;
    %       for i = 2:nnode_ext(1) - 1
    %         i0 = i - 1;
    %
    %         mass = vol_of_node * (mesh.rho_for_2dcell(j0, i0) + mesh.rho_for_2dcell(j0, i) + ...
    %           mesh.rho_for_2dcell(j, i0) + mesh.rho_for_2dcell(j, i)) + tiny;
    %
    %         if direction == 1
    %           obj.force_for_2dnode(j, i) = ...
    %             mesh.pres_for_2dcell(j0, i0) + mesh.pres_for_2dcell(j, i0) - ...
    %             mesh.pres_for_2dcell(j0, i ) - mesh.pres_for_2dcell(j, i) + ...
    %             obj.qvis_for_2dcell(j0, i0) + obj.qvis_for_2dcell(j, i0) - ...
    %             obj.qvis_for_2dcell(j0, i ) - obj.qvis_for_2dcell(j, i);
    %         elseif direction == 2
    %           obj.force_for_2dnode(j, i) = ...
    %             mesh.pres_for_2dcell(j0, i0) + mesh.pres_for_2dcell(j0, i) - ...
    %             mesh.pres_for_2dcell(j , i0) - mesh.pres_for_2dcell(j, i) + ...
    %             obj.qvis_for_2dcell(j0, i0) + obj.qvis_for_2dcell(j0, i) - ...
    %             obj.qvis_for_2dcell(j , i0) - obj.qvis_for_2dcell(j, i);
    %         end
    %
    %         obj.force_for_2dnode(j, i) = obj.force_for_2dnode(j, i) * (area_of_corner / mass);
    %       end
    %     end
    %
    %   elseif dim == 3
    %     % Initialize the force_for_3dnode array if it is empty
    %     if isempty(obj.force_for_3dnode)
    %       obj.force_for_3dnode = zeros(nnode_ext(3), nnode_ext(2), nnode_ext(1));
    %     end
    %
    %     % Calculate force for 3D nodes
    %     for k = 2:nnode_ext(3) - 1
    %       k0 = k - 1;
    %       for j = 2:nnode_ext(2) - 1
    %         j0 = j - 1;
    %         for i = 2:nnode_ext(1) - 1
    %           i0 = i - 1;
    %
    %           mass = vol_of_node * (mesh.rho_for_3dcell(k0, j0, i0) + mesh.rho_for_3dcell(k0, j0, i) + ...
    %             mesh.rho_for_3dcell(k0, j, i0) + mesh.rho_for_3dcell(k, j0, i0) + ...
    %             mesh.rho_for_3dcell(k0, j, i) + mesh.rho_for_3dcell(k, j0, i) + ...
    %             mesh.rho_for_3dcell(k, j, i0) + mesh.rho_for_3dcell(k, j, i)) + tiny;
    %
    %           if direction == 1
    %             obj.force_for_3dnode(k, j, i) = (mesh.pres_for_3dcell(k0, j0, i0) + mesh.pres_for_3dcell(k0, j, i0) + ...
    %               mesh.pres_for_3dcell(k, j0, i0) + mesh.pres_for_3dcell(k, j, i0) - ...
    %               mesh.pres_for_3dcell(k0, j0, i) - mesh.pres_for_3dcell(k0, j, i) - ...
    %               mesh.pres_for_3dcell(k, j0, i) - mesh.pres_for_3dcell(k, j, i)) + ...
    %               (obj.qvis_for_3dcell(k0, j0, i0) + obj.qvis_for_3dcell(k0, j, i0) + ...
    %               obj.qvis_for_3dcell(k, j0, i0) + obj.qvis_for_3dcell(k, j, i0) - ...
    %               obj.qvis_for_3dcell(k0, j0, i) - obj.qvis_for_3dcell(k0, j, i) - ...
    %               obj.qvis_for_3dcell(k, j0, i) - obj.qvis_for_3dcell(k, j, i));
    %           elseif direction == 2
    %             obj.force_for_3dnode(k, j, i) = (mesh.pres_for_3dcell(k0, j0, i0) + mesh.pres_for_3dcell(k0, j0, i) + ...
    %               mesh.pres_for_3dcell(k, j0, i0) + mesh.pres_for_3dcell(k, j0, i) - ...
    %               mesh.pres_for_3dcell(k0, j, i0) - mesh.pres_for_3dcell(k0, j, i) - ...
    %               mesh.pres_for_3dcell(k, j, i0) - mesh.pres_for_3dcell(k, j, i)) + ...
    %               (obj.qvis_for_3dcell(k0, j0, i0) + obj.qvis_for_3dcell(k0, j0, i) + ...
    %               obj.qvis_for_3dcell(k, j0, i0) + obj.qvis_for_3dcell(k, j0, i) - ...
    %               obj.qvis_for_3dcell(k0, j, i0) - obj.qvis_for_3dcell(k0, j, i) - ...
    %               obj.qvis_for_3dcell(k, j, i0) - obj.qvis_for_3dcell(k, j, i));
    %           elseif direction == 3
    %             obj.force_for_3dnode(k, j, i) = (mesh.pres_for_3dcell(k0, j0, i0) + mesh.pres_for_3dcell(k0, j0, i) + ...
    %               mesh.pres_for_3dcell(k0, j, i0) + mesh.pres_for_3dcell(k0, j, i) - ...
    %               mesh.pres_for_3dcell(k, j0, i0) - mesh.pres_for_3dcell(k, j0, i) - ...
    %               mesh.pres_for_3dcell(k, j, i0) - mesh.pres_for_3dcell(k, j, i)) + ...
    %               (obj.qvis_for_3dcell(k0, j0, i0) + obj.qvis_for_3dcell(k0, j0, i) + ...
    %               obj.qvis_for_3dcell(k0, j, i0) + obj.qvis_for_3dcell(k0, j, i) - ...
    %               obj.qvis_for_3dcell(k, j0, i0) - obj.qvis_for_3dcell(k, j0, i) - ...
    %               obj.qvis_for_3dcell(k, j, i0) - obj.qvis_for_3dcell(k, j, i));
    %           end
    %
    %           obj.force_for_3dnode(k, j, i) = obj.force_for_3dnode(k, j, i) * (area_of_corner / mass);
    %         end
    %       end
    %     end
    %   end
    % end
    %% New compute_force function
    function compute_force(obj, mesh, direction)
      % compute_force - Calculate forces on the mesh
      % Inputs:
      %   obj       - Instance of the c_update class
      %   mesh      - Mesh object containing problem properties
      %   direction - Direction of the force calculation (1, 2, or 3 for x, y, z)

      global tiny % Global variable for small values

      % Extract required mesh properties
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      dx = mesh.dx(1:dim);

      % Initialize variables
      nnode_ext = ncell + 2 * nbdry + 1; % Extended nodes in each dimension
      vol_of_corner = prod(0.5 * dx);   % Volume of a corner
      area_of_node = prod(dx) / dx(direction); % Area of a node face perpendicular to the direction

      if dim == 2
        obj.force_for_2dnode = zeros(nnode_ext(2), nnode_ext(1)); % Initialize force array

        if direction == 1
          for j = 2:nnode_ext(2)-1
            j0 = j - 1;
            for i = 2:nnode_ext(1)-1
              i0 = i - 1;

              % Compute mass for the node
              mass = vol_of_corner * (mesh.rho_for_2dcell(j0, i0) + mesh.rho_for_2dcell(j0, i) + ...
                mesh.rho_for_2dcell(j, i0) + mesh.rho_for_2dcell(j, i)) + tiny;

              % Compute force in x-direction
              obj.force_for_2dnode(j, i) = ...
                (mesh.pres_for_2dcell(j0, i0) + mesh.pres_for_2dcell(j, i0) - ...
                mesh.pres_for_2dcell(j0, i) - mesh.pres_for_2dcell(j, i)) + ...
                (obj.qvis_for_2dcell(j0, i0) + obj.qvis_for_2dcell(j, i0) - ...
                obj.qvis_for_2dcell(j0, i) - obj.qvis_for_2dcell(j, i));

              % Apply scaling
              obj.force_for_2dnode(j, i) = 0.5 * obj.force_for_2dnode(j, i) * (area_of_node / mass);
            end
          end
        elseif direction == 2
          for j = 2:nnode_ext(2)-1
            j0 = j - 1;
            for i = 2:nnode_ext(1)-1
              i0 = i - 1;

              % Compute mass for the node
              mass = vol_of_corner * (mesh.rho_for_2dcell(j0, i0) + mesh.rho_for_2dcell(j0, i) + ...
                mesh.rho_for_2dcell(j, i0) + mesh.rho_for_2dcell(j, i)) + tiny;

              % Compute force in y-direction
              obj.force_for_2dnode(j, i) = ...
                (mesh.pres_for_2dcell(j0, i0) + mesh.pres_for_2dcell(j0, i) - ...
                mesh.pres_for_2dcell(j, i0) - mesh.pres_for_2dcell(j, i)) + ...
                (obj.qvis_for_2dcell(j0, i0) + obj.qvis_for_2dcell(j0, i) - ...
                obj.qvis_for_2dcell(j, i0) - obj.qvis_for_2dcell(j, i));

              % Apply scaling
              obj.force_for_2dnode(j, i) = 0.5 * obj.force_for_2dnode(j, i) * (area_of_node / mass);
            end
          end
        end
      elseif dim == 3
        obj.force_for_3dnode = zeros(nnode_ext(3), nnode_ext(2), nnode_ext(1)); % Initialize force array


        for k = 2:nnode_ext(3)-1
          k0 = k - 1;
          for j = 2:nnode_ext(2)-1
            j0 = j - 1;
            for i = 2:nnode_ext(1)-1
              i0 = i - 1;

              % Compute mass for the node
              mass = vol_of_corner * (mesh.rho_for_3dcell(k0, j0, i0) + mesh.rho_for_3dcell(k0, j0, i) + ...
                mesh.rho_for_3dcell(k0, j, i0) + mesh.rho_for_3dcell(k, j0, i0) + ...
                mesh.rho_for_3dcell(k0, j, i) + mesh.rho_for_3dcell(k, j0, i) + ...
                mesh.rho_for_3dcell(k, j, i0) + mesh.rho_for_3dcell(k, j, i)) + tiny;

              % Compute force for the specified direction
              if direction == 1
                obj.force_for_3dnode(k, j, i) = ...
                  (mesh.pres_for_3dcell(k0, j0, i0) + mesh.pres_for_3dcell(k0, j, i0) + ...
                  mesh.pres_for_3dcell(k, j0, i0) + mesh.pres_for_3dcell(k, j, i0) - ...
                  mesh.pres_for_3dcell(k0, j0, i) - mesh.pres_for_3dcell(k0, j, i) - ...
                  mesh.pres_for_3dcell(k, j0, i) - mesh.pres_for_3dcell(k, j, i)) + ...
                  (obj.qvis_for_3dcell(k0, j0, i0) + obj.qvis_for_3dcell(k0, j, i0) + ...
                  obj.qvis_for_3dcell(k, j0, i0) + obj.qvis_for_3dcell(k, j, i0) - ...
                  obj.qvis_for_3dcell(k0, j0, i) - obj.qvis_for_3dcell(k0, j, i) - ...
                  obj.qvis_for_3dcell(k, j0, i) - obj.qvis_for_3dcell(k, j, i));
              elseif direction == 2
                obj.force_for_3dnode(k, j, i) = ...
                  (mesh.pres_for_3dcell(k0, j0, i0) + mesh.pres_for_3dcell(k0, j0, i) + ...
                  mesh.pres_for_3dcell(k, j0, i0) + mesh.pres_for_3dcell(k, j0, i) - ...
                  mesh.pres_for_3dcell(k0, j, i0) - mesh.pres_for_3dcell(k0, j, i) - ...
                  mesh.pres_for_3dcell(k, j, i0) - mesh.pres_for_3dcell(k, j, i)) + ...
                  (obj.qvis_for_3dcell(k0, j0, i0) + obj.qvis_for_3dcell(k0, j0, i) + ...
                  obj.qvis_for_3dcell(k, j0, i0) + obj.qvis_for_3dcell(k, j0, i) - ...
                  obj.qvis_for_3dcell(k0, j, i0) - obj.qvis_for_3dcell(k0, j, i) - ...
                  obj.qvis_for_3dcell(k, j, i0) - obj.qvis_for_3dcell(k, j, i));
              elseif direction == 3
                obj.force_for_3dnode(k, j, i) = ...
                  (mesh.pres_for_3dcell(k0, j0, i0) + mesh.pres_for_3dcell(k0, j0, i) + ...
                  mesh.pres_for_3dcell(k0, j, i0) + mesh.pres_for_3dcell(k0, j, i) - ...
                  mesh.pres_for_3dcell(k, j0, i0) - mesh.pres_for_3dcell(k, j0, i) - ...
                  mesh.pres_for_3dcell(k, j, i0) - mesh.pres_for_3dcell(k, j, i)) + ...
                  (obj.qvis_for_3dcell(k0, j0, i0) + obj.qvis_for_3dcell(k0, j0, i) + ...
                  obj.qvis_for_3dcell(k0, j, i0) + obj.qvis_for_3dcell(k0, j, i) - ...
                  obj.qvis_for_3dcell(k, j0, i0) - obj.qvis_for_3dcell(k, j0, i) - ...
                  obj.qvis_for_3dcell(k, j, i0) - obj.qvis_for_3dcell(k, j, i));
              end

              % Apply scaling
              obj.force_for_3dnode(k, j, i) = 0.25 * obj.force_for_3dnode(k, j, i) * (area_of_node / mass);
            end
          end
        end
      end
    end


    %%
    function update_vel_comp(obj, mesh, dt, dir)
      % Update the velocity components for the mesh
      % Inputs:
      %   - obj: instance of c_update
      %   - mesh: instance of the mesh class with properties defining the problem
      %   - dt: time step size
      %   - dir: direction of velocity update (0, 1, or 2)

      % Some shortcuts
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      hdt = 0.5 * dt;  % Half of the time step

      nnode_ext = ncell + 2 * nbdry + 1;  % Extended number of nodes

      % Check dimensionality and update velocity components
      if dim == 2
        % Update velocity for 2D nodes
        for j = 2:nnode_ext(2)-1
          for i = 2:nnode_ext(1)-1
            if i == 53
              disp("");
            end
            hdv = hdt * obj.force_for_2dnode(j, i);
            mesh.vav_for_2dnode(j, i, dir) = mesh.vel_for_2dnode(j, i, dir) + hdv;
            mesh.vel_for_2dnode(j, i, dir) = mesh.vel_for_2dnode(j, i, dir) + 2 * hdv;
          end
        end

      elseif dim == 3
        % Update velocity for 3D nodes
        for k = 2:nnode_ext(3)-1
          for j = 2:nnode_ext(2)-1
            for i = 2:nnode_ext(1)-1
              hdv = hdt * obj.force_for_3dnode(k, j, i);
              mesh.vav_for_3dnode(k, j, i, dir) = mesh.vel_for_3dnode(k, j, i, dir) + hdv;
              mesh.vel_for_3dnode(k, j, i, dir) = mesh.vel_for_3dnode(k, j, i, dir) + 2 * hdv;
            end
          end
        end
      end
    end

    function courant = update_density(obj, mesh, dt)
      % Update the density in the mesh and compute the Courant number
      % Inputs:
      %   - obj: instance of c_update
      %   - mesh: instance of the mesh class with properties defining the problem
      %   - mat: instance of the mat class with material-related properties
      %   - dt: time step size
      % Output:
      %   - cournt: maximum Courant number

      % Some shorthands
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      nmat = mesh.nmat_mesh;
      courant = 0.0;
      ncell_ext = ncell + 2 * nbdry;

      if dim == 2
        % Loop through 2D cells, skipping boundaries
        for j = 2:ncell_ext(2)-1
          for i = 2:ncell_ext(1)-1
            % Reset the cell density
            mesh.rho_for_2dcell(j, i) = 0.0;

            % Update Courant number
            courant = max(courant, -obj.divu_for_2dcell(j, i) * dt);

            % Compute factor for updating density
            factor = 1.0 / (1.0 + obj.divu_for_2dcell(j, i) * dt);

            % Update material densities and cell density
            for m = 1:nmat
              mesh.rho_2dmat(j, i, m) = mesh.rho_2dmat(j, i, m) * factor;
              mesh.rho_for_2dcell(j, i) = mesh.rho_for_2dcell(j, i) ...
                + (mesh.vf_2dmat(j, i, m) *mesh. rho_2dmat(j, i, m));
            end
          end
        end

      elseif dim == 3
        % Loop through 3D cells, skipping boundaries
        for k = 2:ncell_ext(3)-1
          for j = 2:ncell_ext(2)-1
            for i = 2:ncell_ext(1)-1
              % Reset the cell density
              mesh.rho_for_3dcell(k, j, i) = 0.0;

              % Update Courant number
              courant = max(courant, -obj.divu_for_3dcell(k, j, i) * dt);

              % Compute factor for updating density
              factor = 1.0 / (1.0 + obj.divu_for_3dcell(k, j, i) * dt);

              % Update material densities and cell density
              for m = 1:nmat
                mesh.rho_3dmat(k, j, i, m) = mesh.rho_3dmat(k, j, i, m) * factor;
                mesh.rho_for_3dcell(k, j, i) = mesh.rho_for_3dcell(k, j, i) ...
                  + (mesh.vf_3dmat(k, j, i, m) * mesh.rho_3dmat(k, j, i, m));
              end
            end
          end
        end
      end

    end

    function update_energy(obj, mesh, dt)
      % Update the internal energy in the mesh
      % Inputs:
      %   - obj: instance of c_update
      %   - mesh: instance of the mesh class with properties defining the problem
      %   - mat: instance of the mat class with material-related properties
      %   - dt: time step size
      %   - rho_in_mixcell_old, es_in_mixcell_old: previous state variables for mixed cells
      global tiny;

      % Some shorthands
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      ncell_ext = ncell + 2 * nbdry;
      nmat = mesh.nmat_mesh;
      % Handling the 2D case
      if dim == 2
        % Loop through each cell and material for 2D
        for j = 1:ncell_ext(2)
          for i = 1:ncell_ext(1)
            for m = 1:nmat
              mesh.ei_2dmat(j,i,m) = mesh.ei_2dmat(j,i,m) / (mesh.rho_2dmat(j,i,m) + tiny); % specific internal energy
            end
          end
        end

        % Adjust the specific internal energy for 2D cells considering pressure and divergence
        for j = 2:ncell_ext(2)-1
          for i = 2:ncell_ext(1)-1
            for m = 1:nmat
              work = (mesh.pres_2dmat(j,i,m) / (mesh.rho_for_2dcell(j,i) + tiny)) ...
                * obj.divu_for_2dcell(j,i) * dt;
              mesh.ei_2dmat(j,i,m) = mesh.ei_2dmat(j,i,m) - work;  % Update specific internal energy
            end
          end
        end

        % Handling the 3D case
      elseif dim == 3
        % Loop through each cell and material for 3D
        for k = 1:ncell_ext(3)
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              for m = 1:nmat
                mesh.ei_3dmat(k,j,i,m) = mesh.ei_3dmat(k,j,i,m) / (mesh.rho_3dmat(k,j,i,m) + tiny); % specific internal energy
              end
            end
          end
        end

        % Adjust the specific internal energy for 3D cells considering pressure and divergence
        for k = 2:ncell_ext(3)-1
          for j = 2:ncell_ext(2)-1
            for i = 2:ncell_ext(1)-1
              for m = 1:nmat
                work = (mesh.pres_3dmat(k,j,i,m) / (mesh.rho_for_3dcell(k,j,i) + tiny)) ...
                  * obj.divu_for_3dcell(k,j,i) * dt;
                mesh.ei_3dmat(k,j,i,m) = mesh.ei_3dmat(k,j,i,m) - work;  % Update specific internal energy
              end
            end
          end
        end
      end
    end

    %% pressure
    function update_pressure(~, mesh, is_solid_ea_mat, gamma_ea_mat)
      % update_pressure Updates pressure and internal energy for each cell in mesh.
      %
      % Inputs:
      %   - mesh: Struct containing 2D or 3D mesh data.
      %   - nmat: Number of materials.
      %   - is_solid_ea_mat: Logical array indicating solid materials.
      %   - gamma_ea_mat: Array of gamma values for each material.
      global eos;

      dim = mesh.dim_prob;
      nbdry = mesh.nbdry_prob;
      ncell_ext = mesh.ncell_prob + 2 * nbdry;
      nmat = mesh.nmat_mesh;
      if dim == 2
        % Update for 2D cells
        for j = 1:ncell_ext(2)
          for i = 1:ncell_ext(1)
            mesh.ei_for_2dcell(j, i) = 0.0;
            mesh.pres_for_2dcell(j, i) = 0.0;
            for m = 1:nmat
              % Update specific internal energy density
              mesh.ei_2dmat(j, i, m) = mesh.ei_2dmat(j, i, m) * mesh.rho_2dmat(j, i, m);

              % Update cell-level energy and pressure
              mesh.ei_for_2dcell(j, i) = mesh.ei_for_2dcell(j, i) + ...
                mesh.vf_2dmat(j, i, m) * mesh.ei_2dmat(j, i, m);

              if ~is_solid_ea_mat(m)
                % Calculate pressure for non-solid materials
                mesh.pres_2dmat(j, i, m) = (gamma_ea_mat(m) - 1.0) * mesh.ei_2dmat(j, i, m);
              else
                % Use Mie-Gruneisen EOS for solid materials
                mesh.pres_2dmat(j, i, m) = eos.p_mie_gruneisen(mesh.rho_2dmat(j, i, m), ...
                  mesh.ei_2dmat(j, i, m));
              end

              mesh.pres_for_2dcell(j, i) = mesh.pres_for_2dcell(j, i) + ...
                mesh.vf_2dmat(j, i, m) * mesh.pres_2dmat(j, i, m);
            end
          end
        end
      elseif dim == 3
        % Update for 3D cells
        for k = 1:ncell_ext(3)
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              mesh.ei_for_3dcell(k, j, i) = 0.0;
              mesh.pres_for_3dcell(k, j, i) = 0.0;
              for m = 1:nmat
                % Update specific internal energy density
                mesh.ei_3dmat(k, j, i, m) = mesh.ei_3dmat(k, j, i, m) * mesh.rho_3dmat(k, j, i, m);

                % Update cell-level energy and pressure
                mesh.ei_for_3dcell(k, j, i) = mesh.ei_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.ei_3dmat(k, j, i, m);

                if ~is_solid_ea_mat(m)
                  % Calculate pressure for non-solid materials
                  mesh.pres_3dmat(k, j, i, m) = (gamma_ea_mat(m) - 1.0) * mesh.ei_3dmat(k, j, i, m);
                else
                  % Use Mie-Gruneisen EOS for solid materials
                  mesh.pres_3dmat(k, j, i, m) = eos.p_mie_gruneisen(mesh.rho_3dmat(k, j, i, m), ...
                    mesh.ei_3dmat(k, j, i, m));
                end

                mesh.pres_for_3dcell(k, j, i) = mesh.pres_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.pres_3dmat(k, j, i, m);
              end
            end
          end
        end
      end
    end

  end
end