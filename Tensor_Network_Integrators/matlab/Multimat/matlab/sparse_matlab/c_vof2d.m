classdef c_vof2d < handle
  % Volume fraction in 2D

  properties (Constant)
    accuracy = 1.0e-10;
    niter_mx = 50;
  end
  methods

    %%%%%%%%%%%%%%% implemented but not tested
    function [nnode_new, coords_new, nnode_interface, nodes_interface, ...
        nnode_lower, nodelist_lower, nnode_upper, nodelist_upper] = ...
        find_interface2d(~, nnode, coords, nodelist, norm, distance, node_loc)
      % find_interface2d - Find the intersection of a 2D interface with a polygon.
      %
      % This function determines the intersection of a 2D interface with a given
      % polygon. It classifies nodes as being on, above, or below the interface,
      % and updates the corresponding node lists and coordinates based on this
      % classification.
      %
      % Inputs:
      %   nnode            - (int) Number of nodes in the polygon.
      %   coords           - (double array) Coordinates of the polygon nodes,
      %                      size: [2*nnode, 1].
      %   nodelist         - (int array) List of node indices for the polygon,
      %                      size: [nnode, 1].
      %   norm             - (double array) Normal vector defining the interface,
      %                      size: [2, 1].
      %   distance         - (double) Distance of the interface from the origin.
      %   node_loc         - (int array) Classification of nodes relative to the
      %                      interface.
      %                      Values: 0 (on interface), -1 (below), 1 (above),
      %                      size: [nnode, 1].
      %
      % Outputs:
      %   nnode_new        - (int) Number of new nodes created at intersections.
      %   coords_new       - (double array) Coordinates of new nodes,
      %                      size: [2*(nnode_new), 1].
      %   nnode_interface  - (int) Number of nodes on the interface.
      %   nodes_interface  - (int array) Indices of nodes on the interface,
      %                      size: [2, 1].
      %   nnode_lower      - (int) Number of nodes below the interface.
      %   nodelist_lower   - (int array) List of node indices below the interface.
      %   nnode_upper      - (int) Number of nodes above the interface.
      %   nodelist_upper   - (int array) List of node indices above the interface.
      %
      %
      % Original C declaration:
      % void find_interface2d(int nnode, double *coords, int *nodelist,
      %                       double *norm, double distance,
      %                       int *node_loc,
      %                       int *nnode_new, double *coords_new,
      %                       int *nnode_interface, int *nodes_interface,
      %                       int *nnode_lower, int **nodelist_lower,
      %                       int *nnode_upper, int **nodelist_upper)

      % Initialize outputs
      nnode_new = 0;
      nnode_interface = 0;
      nnode_lower = 0;
      nnode_upper = 0;

      % Reshape coords
      if size(coords, 2) == 2
        coords = coords';
        coords = coords(:);
      end

      % Preallocate arrays for new nodes and interface nodes
      coords_new = zeros(4, 1);
      nodes_interface = zeros(4, 1);
      nodelist_lower = [];
      nodelist_upper = [];

      % Allocate working arrays
      nodelist_low = zeros(nnode, 1);
      nodelist_up = zeros(nnode, 1);

      for i0 = 1:nnode
        i1 = mod(i0, nnode) + 1;
        n1 = nodelist(i1);
        n0 = nodelist(i0);

        if node_loc(n0) == 0 && node_loc(n1) == 0
          % Both points are on the line
          nodelist_low(nnode_lower + 1) = n0;
          nodelist_low(nnode_lower + 2) = n1;
          nnode_lower = nnode_lower + 2;
          nodelist_up(nnode_upper + 1) = n0;
          nodelist_up(nnode_upper + 2) = n1;
          nnode_upper = nnode_upper + 2;

          nodes_interface(1) = n0;
          nodes_interface(2) = n1;
          nnode_interface = 2;

        elseif node_loc(n0) == 0
          nodelist_low(nnode_lower + 1) = n0;
          nodelist_up(nnode_upper + 1) = n0;
          nnode_lower = nnode_lower + 1;
          nnode_upper = nnode_upper + 1;

          if nnode_interface > 0
            if nodes_interface(nnode_interface) ~= n0
              nnode_interface = nnode_interface + 1;
              nodes_interface(nnode_interface) = n0;
            end
          else
            nodes_interface(1) = n0;
            nnode_interface = 1;
          end

          if node_loc(n1) < 0
            nodelist_low(nnode_lower + 1) = n1;
            nnode_lower = nnode_lower + 1;
          else
            nodelist_up(nnode_upper + 1) = n1;
            nnode_upper = nnode_upper + 1;
          end

        elseif node_loc(n1) == 0
          if node_loc(n0) < 0
            nodelist_low(nnode_lower + 1) = n0;
            nnode_lower = nnode_lower + 1;
          elseif node_loc(n0) > 0
            nodelist_up(nnode_upper + 1) = n0;
            nnode_upper = nnode_upper + 1;
          end
          nodelist_low(nnode_lower + 1) = n1;
          nodelist_up(nnode_upper + 1) = n1;
          nnode_lower = nnode_lower + 1;
          nnode_upper = nnode_upper + 1;

          if nnode_interface < 2
            if nnode_interface > 0
              if nodes_interface(nnode_interface) ~= n1
                nnode_interface = nnode_interface + 1;
                nodes_interface(nnode_interface) = n1;
              end
            else
              nodes_interface(1) = n1;
              nnode_interface = 1;
            end
          end

        elseif node_loc(n0) > 0 && node_loc(n1) > 0
          nodelist_up(nnode_upper + 1) = n0;
          nodelist_up(nnode_upper + 2) = n1;
          nnode_upper = nnode_upper + 2;

        elseif node_loc(n0) < 0 && node_loc(n1) < 0
          nodelist_low(nnode_lower + 1) = n0;
          nodelist_low(nnode_lower + 2) = n1;
          nnode_lower = nnode_lower + 2;

        else
          % Find the intersection
          c0 = coords(2*n0-1:2*n0);
          c1 = coords(2*n1-1:2*n1);
          t = (distance - (norm(1) * c0(1) + norm(2) * c0(2))) / ...
            (norm(1) * (c1(1) - c0(1)) + norm(2) * (c1(2) - c0(2)));
          t = max(0.0, min(1.0, t));
          nnode_new = nnode_new + 1;
          coords_new(2*nnode_new - 1) = c0(1) + (c1(1) - c0(1))*t;
          coords_new(2*nnode_new)     = c0(2) + (c1(2) - c0(2))*t;

          if node_loc(n0) < 0
            nodelist_low(nnode_lower + 1) = n0;
            nodelist_low(nnode_lower + 2) = nnode + nnode_new;
            nnode_lower = nnode_lower + 2;

            nodelist_up(nnode_upper + 1) = nnode + nnode_new;
            nodelist_up(nnode_upper + 2) = n1;
            nnode_upper = nnode_upper + 2;
          else % if node_loc(n1) < 0
            nodelist_low(nnode_lower + 1) = nnode + nnode_new;
            nodelist_low(nnode_lower + 2) = n1;
            nnode_lower = nnode_lower + 2;

            nodelist_up(nnode_upper + 1) = n0;
            nodelist_up(nnode_upper + 2) = nnode + nnode_new;
            nnode_upper = nnode_upper + 2;
          end

          nnode_interface = nnode_interface + 1;
          nodes_interface(nnode_interface) = nnode + nnode_new;
        end
      end

      % Adjust lower and upper nodelists
      if nnode_lower > 0
        if nodelist_low(1) == nodelist_low(nnode_lower)
          nnode_lower = nnode_lower - 1;
        end
      end

      if nnode_upper > 0
        if nodelist_up(1) == nodelist_up(nnode_upper)
          nnode_upper = nnode_upper - 1;
        end
      end

      % Remove redundant nodes in lower list
      nn = 2;
      for i0 = 2:nnode_lower
        if nodelist_low(i0) ~= nodelist_low(nn - 1)
          nodelist_low(nn) = nodelist_low(i0);
          nn = nn + 1;
        end
      end
      nnode_lower = nn - 1;
      if nnode_lower < 3
        nnode_lower = 0;
        nodelist_lower = [];
      else
        nodelist_lower = nodelist_low(1:nnode_lower);
      end

      % Remove redundant nodes in upper list
      nn = 2;
      for i0 = 2:nnode_upper
        if nodelist_up(i0) ~= nodelist_up(nn - 1)
          nodelist_up(nn) = nodelist_up(i0);
          nn = nn + 1;
        end
      end
      nnode_upper = nn - 1;
      if nnode_upper < 3
        nnode_upper = 0;
        nodelist_upper = [];
      else
        nodelist_upper = nodelist_up(1:nnode_upper);
      end

      % Trim interface list
      nodes_interface = nodes_interface(1:nnode_interface);
      coords_new = coords_new(1:2*nnode_new);
    end

    function [ds_lower, ds_upper, vf_lower, vf_upper, nnode_new, coords_new, ...
        vol_matched, nnode_interface, nodes_interface, nnode_lower, ...
        nodelist_lower, nnode_upper, nodelist_upper] = ...
        bounds_2d(obj, geop, vf_to_match, volume, norm, nnode, ...
        coords, nodelist, node_order_for_ds, ds_ea_node)
      % bounds_2d - Determine the bounds for a 2D volume fraction interface within a polygon.
      %
      % This function calculates the bounds (lower and upper distances) along a specified
      % normal vector for a 2D polygon. It divides the polygon into regions above, below,
      % and on a volume fraction interface and determines which portion of the polygon
      % corresponds to a target volume fraction (`vf_to_match`). The interface nodes and
      % coordinates, along with the corresponding volume fractions, are returned.
      %
      % Inputs:
      %   geop             - (int) Type of geometry to use: 1 for polygonal,
      %                      2 for rz-plane.
      %   vf_to_match      - (double) The target volume fraction to match.
      %   volume           - (double) The total volume of the domain.
      %   norm             - (double array) Normal vector defining the interface,
      %                      size: [2, 1].
      %   nnode            - (int) The number of nodes defining the polygon.
      %   coords           - (double array) Coordinates of the nodes, provided as
      %                      a flattened array where each point has two values (x, y),
      %                      size: [2*nnode, 1].
      %   nodelist         - (int array) List of node indices for the polygon,
      %                      size: [nnode, 1].
      %   node_order_for_ds - (int array) Node order based on distance, size: [nnode, 1].
      %   ds_ea_node       - (double array) Distance for each node, size: [nnode, 1].
      %
      % Outputs:
      %   ds_lower         - (double) Lower bound of the distance along the normal.
      %   ds_upper         - (double) Upper bound of the distance along the normal.
      %   vf_lower         - (double) Volume fraction corresponding to `ds_lower`.
      %   vf_upper         - (double) Volume fraction corresponding to `ds_upper`.
      %   nnode_new        - (int) Number of new nodes created at the interface.
      %   coords_new       - (double array) Coordinates of new nodes created at the
      %                      interface, size: [2*nnode_new, 1].
      %   vol_matched      - (int) Flag indicating if the volume fraction was matched.
      %                      (1 if matched, 0 otherwise).
      %   nnode_interface  - (int) Number of nodes on the interface.
      %   nodes_interface  - (int array) Indices of nodes on the interface, size: [2, 1].
      %   nnode_lower      - (int) Number of nodes below the interface.
      %   nodelist_lower   - (int array) List of node indices below the interface,
      %                      size: [nnode_lower, 1].
      %   nnode_upper      - (int) Number of nodes above the interface.
      %   nodelist_upper   - (int array) List of node indices above the interface,
      %                      size: [nnode_upper, 1].
      %
      % Function Logic:
      %   1. The function starts by ordering the nodes of the polygon based on their distance
      %      along the normal vector, grouping nodes with equal distances.
      %   2. It then iterates through these groups of nodes, adjusting the interface distance
      %      (`ds`) to find the interface position that matches the target volume fraction
      %      (`vf_to_match`).
      %   3. For each potential interface, it divides the nodes into three categories:
      %      - Nodes below the interface,
      %      - Nodes above the interface,
      %      - Nodes on the interface.
      %   4. The polygon is partitioned into these regions and the volume fraction is calculated
      %      based on the lower region.
      %   5. The algorithm performs a search, adjusting `ds` (the interface position) to match
      %      the volume fraction (`vf_to_match`) with a specified accuracy.
      %   6. If the volume fraction is matched within the desired tolerance, the function sets the
      %      lower and upper bounds (`ds_lower`, `ds_upper`), and returns the interface nodes.
      %   7. If the volume fraction cannot be exactly matched, it returns the bounds that bracket
      %      the closest approximation of the volume fraction.
      %
      % Notes:
      %   - The function uses `find_interface2d` to calculate the intersection between the
      %     interface and the polygon.
      %   - This function assumes that the `coords_new` and `nodes_interface` arrays have been
      %     preallocated.
      %
      %
      % Original C declaration:
      % void bounds_2d(int geop,
      %                double vf_to_match, double volume, double *norm,
      %                int nnode, double *coords, int *nodelist,
      %                int *node_order_for_ds, double *ds_ea_node,
      %                double *ds_lower, double *ds_upper,
      %                double *vf_lower, double *vf_upper,
      %                int *nnode_new, double *coords_new,
      %                int *vol_matched,
      %                int *nnode_interface, int *nodes_interface,
      %                int *nnode_lower, int **nodelist_lower,
      %                int *nnode_upper, int **nodelist_upper)

      % Globals
      global util;

      % Initialize outputs
      ds_lower = 0.0;
      ds_upper = 0.0;
      vf_lower = 0.0;
      vf_upper = 0.0;
      nnode_new = 0;
      coords_new = [];
      vol_matched = false;
      nnode_interface = 0;
      nodes_interface = zeros(2, 1);
      nnode_lower = 0;
      nodelist_lower = [];
      nnode_upper = 0;
      nodelist_upper = [];

      % Initialize local variables
      dim = 2;
      szdim = dim * 8;  % Size of double in bytes
      ds_ordered = zeros(1, nnode);
      node_loc = zeros(1, nnode);
      nnode_ea_ds = zeros(1, nnode);
      nodelist_ea_ds = zeros(1, nnode);
      coords_work = zeros(2, nnode);

      % Initial setup
      n = node_order_for_ds(1);
      ds = ds_ea_node(n);
      ds_ordered(1) = ds;
      nnode_ea_ds(1) = 1;
      nodelist_ea_ds(1) = n;
      offset = 1;
      ndistance = 1;
      idx_next = 2;
      found = true;

      while found
        found = false;
        for i = idx_next:nnode
          n = node_order_for_ds(i);
          if ds_ea_node(n) - ds > obj.accuracy
            ds = ds_ea_node(n);
            ds_ordered(ndistance + 1) = ds;
            nnode_ea_ds(ndistance + 1) = 1;
            nodelist_ea_ds(offset + 1) = n;
            ndistance = ndistance + 1;
            offset = offset + 1;
            idx_next = i + 1;
            found = true;
            break;
          else
            nnode_ea_ds(ndistance) = nnode_ea_ds(ndistance) + 1;
            nodelist_ea_ds(offset + 1) = n;
            offset = offset + 1;
          end
        end
      end

      idx = 1;
      ds_lower = ds_ordered(1);
      vf_lower = 0.0;

      while idx < ndistance
        ds = ds_ordered(idx + 1);

        % Reset nodelist_lower and nodelist_upper
        nodelist_lower = [];
        nodelist_upper = [];

        % Mark all nodes as initially above the interface
        node_loc(:) = 1;

        % Mark nodes below the interface
        offset = sum(nnode_ea_ds(1:idx));
        node_loc(nodelist_ea_ds(1:offset)) = -1;

        % Mark nodes on the interface
        nn = nnode_ea_ds(idx + 1);
        nodes_interface(1:nn) = nodelist_ea_ds(offset + 1:offset + nn);
        node_loc(nodes_interface(1:nn)) = 0;

        if idx == ndistance - 1
          % Last distance point
          nnode_new = 0;
          coords_new = [];
          nnode_interface = nn;
          nodes_interface = nodes_interface(1:nn);
          nnode_upper = 0;
          nnode_lower = nnode;
          nodelist_lower = nodelist;
          vf = 1.0;
        else
          % Call find_interface2d method
          [nnode_new, coords_new, nnode_interface, nodes_interface, ...
            nnode_lower, nodelist_lower, nnode_upper, nodelist_upper] = ...
            obj.find_interface2d(nnode, coords, nodelist, norm, ds, node_loc);

          % Calculate the volume of the lower polygon
          vol = 0.0;
          for i = 1:nnode_lower
            n0 = nodelist_lower(i);
            if n0 <= nnode
              c0 = coords(2*n0-1:2*n0);
            else
              n1 = n0 - nnode;
              c0 = coords_new(2*n1-1:2*n1);
            end
            coords_work(:, i) = c0;
          end

          if geop == 1
            vol = util.cal_poly_area(nnode_lower, coords_work(:), nnode_lower, []);
          elseif geop == 2
            vol = util.rz_area(nnode_lower, coords_work(:));
          end

          vf = vol / volume;
        end

        % Check if the volume fraction matches
        if abs(vf - vf_to_match) < obj.accuracy
          vol_matched = true;
          ds_lower = ds;
          ds_upper = ds;
          vf_lower = vf;
          vf_upper = vf;
          break;
        elseif idx == ndistance - 1
          ds_upper = ds;
          vf_upper = vf;
          break;
        elseif vf > vf_to_match
          ds_upper = ds;
          vf_upper = vf;
          break;
        else
          vf_lower = vf;
          ds_lower = ds;
          idx = idx + 1;
        end
      end
    end

    function [distance, nnode_new, coords_new, nnode_interface, nodes_interface, ...
        nnode_lower, nodelist_lower, nnode_upper, nodelist_upper] = ...
        cal_distance2d(obj, geop, vf_to_match, volume, norm, ...
        nnode, coords, nodelist)
      % cal_distance2d - Calculate the distance along a normal vector for a given
      %                  volume fraction.
      %
      % This function calculates the distance along a normal vector at which the
      % volume fraction of a polygonal region matches a specified target value. The
      % function iteratively refines the distance to achieve the desired volume
      % fraction within a specified accuracy.
      %
      % Inputs:
      %   geop             - (int) Type of geometry to use: 1 for polygonal,
      %                      2 for rz-plane.
      %   vf_to_match      - (double) The target volume fraction to match.
      %   volume           - (double) The total volume of the domain.
      %   norm             - (double array) Normal vector defining the interface,
      %                      size: [2, 1].
      %   nnode            - (int) The number of nodes defining the polygon.
      %   coords           - (double array) Coordinates of the nodes, provided as
      %                      a flattened array where each point has two values (x, y),
      %                      size: [2*nnode, 1].
      %   nodelist         - (int array) List of node indices for the polygon,
      %                      size: [nnode, 1].
      %
      % Outputs:
      %   distance         - (double) The calculated distance along the normal vector
      %                      where the volume fraction matches `vf_to_match`.
      %   nnode_new        - (int) The number of new nodes created at the interface.
      %   coords_new       - (double array) The coordinates of new nodes created at
      %                      the interface, size: [2*nnode_new, 1].
      %   nnode_interface  - (int) The number of nodes on the interface.
      %   nodes_interface  - (int array) The indices of nodes on the interface,
      %                      size: [2, 1].
      %   nnode_lower      - (int) The number of nodes below the interface.
      %   nodelist_lower   - (int array) The list of node indices below the interface,
      %                      size: [nnode_lower, 1].
      %   nnode_upper      - (int) The number of nodes above the interface.
      %   nodelist_upper   - (int array) The list of node indices above the interface,
      %                      size: [nnode_upper, 1].
      %
      % Original C declaration:
      % void cal_distance2d(int geop,
      %                     double vf_to_match, double volume, double *norm,
      %                     int nnode, double *coords, int *nodelist,
      %                     int *nnode_new, double *coords_new,
      %                     double *distance,
      %                     int *nnode_interface, int *nodes_interface,
      %                     int *nnode_lower, int **nodelist_lower,
      %                     int *nnode_upper, int **nodelist_upper)
      global util;

      % Initialize variables
      dim = 2;

      % Allocate working arrays
      coords_work = zeros(2, nnode);
      node_loc = zeros(nnode, 1);

      % Order nodes along the normal vector
      [node_order_for_ds, ds_ea_node] = ...
        util.order_nodes_along_norm(dim, norm, nnode, coords);

      % Calculate bounds
      [ds_lower, ds_upper, vf_lower, vf_upper, nnode_new, coords_new, ...
        vol_matched, nnode_interface, nodes_interface, nnode_lower, ...
        nodelist_lower, nnode_upper, nodelist_upper] = ...
        obj.bounds_2d(geop, vf_to_match, volume, norm, ...
        nnode, coords, nodelist, ...
        node_order_for_ds, ds_ea_node);

      if vol_matched
        distance = ds_upper;
      else
        err = 1.0;
        ds1 = ds_upper;
        ds0 = ds_lower;
        vf0 = vf_lower;
        vf1 = vf_upper;
        vf = vf0;

        distance = 0.5 * (ds0 + ds1);
        niter = 0;
        done = false;

        while ~done && (niter < obj.niter_mx)

          % Reset nodelists
          nodelist_lower = [];
          nodelist_upper = [];

          % Update node locations
          for i = 1:nnode
            if abs(ds_ea_node(i) - distance) <= obj.accuracy
              node_loc(i) = 0;  % On the interface
            elseif ds_ea_node(i) > distance
              node_loc(i) = 1;  % Above the interface
            else
              node_loc(i) = -1; % Below the interface
            end
          end

          % Call find_interface2d method
          [nnode_new, coords_new, nnode_interface, nodes_interface, ...
            nnode_lower, nodelist_lower, nnode_upper, nodelist_upper] = ...
            obj.find_interface2d(nnode, coords, nodelist, norm, ...
            distance, node_loc);

          if nnode_interface < 2
            % Handle case when no valid interface is found
            nnode_interface = 0;
            nodes_interface = [];
            nnode_new = 0;
            coords_new = [];

            if vf_to_match > 0.5
              % All nodes are above the interface
              nnode_lower = 0;
              nodelist_lower = [];
              nnode_upper = nnode;
              nodelist_upper = nodelist;
              n = node_order_for_ds(end);
              distance = ds_ea_node(n);
            else
              % All nodes are below the interface
              nnode_upper = 0;
              nnode_lower = nnode;
              nodelist_upper = [];
              nodelist_lower = nodelist;
              n = node_order_for_ds(1);
              distance = ds_ea_node(n);
            end
            done = true;
          else
            % Calculate the volume of the lower polygon
            for i = 1:nnode_lower
              n = nodelist_lower(i);
              if n <= nnode
                c0 = coords(2*n-1:2*n);
              else
                n = n - nnode;
                c0 = coords_new(2*n-1:2*n);
              end
              coords_work(:, i) = c0;
            end

            if geop == 1
              vol = util.cal_poly_area(nnode_lower, coords_work(:), nnode_lower, []);
            elseif geop == 2
              vol = util.rz_area(nnode_lower, coords_work(:));
            end

            vf = vol / volume;
            err = abs(vf - vf_to_match);

            if (err <= obj.accuracy) || (ds1 - ds0 <= obj.accuracy)
              done = true;
            else
              if vf > vf_to_match
                vf1 = vf;
                ds1 = distance;
              else
                vf0 = vf;
                ds0 = distance;
              end
              distance = 0.5 * (ds0 + ds1);
              niter = niter + 1;
            end
          end
        end

        assert(niter < obj.niter_mx, 'Maximum number of iterations exceeded');
      end
    end

    function [coords_final, nnode_final, nodes_for_interface, ...
        nodelist_for_mpoly, nnode_for_interface, nnode_for_mpoly] = ...
        reconstruct2d_nmat_pagosa(obj, geop, xl, dx, ...
        nmat_mesh, matid_mesh, vf_mesh)
      % reconstruct2d_nmat_pagosa - Reconstruct the 2D geometry for multiple
      %                             materials in a PAGOSA-style mesh.
      %
      % This function reconstructs the 2D geometry of a cell that contains multiple
      % materials using the PAGOSA method. It calculates the material interfaces and
      % reconstructs the polygons representing each material within the cell.
      %
      % Inputs:
      %   geop                - (int) Type of geometry to use: 1 for polygonal,
      %                         2 for rz-plane.
      %   xl                  - (double array) The lower bounds of the cell in the
      %                         x and y directions, size: [2, 1].
      %   dx                  - (double array) The grid spacing in the x and y
      %                         directions, size: [2, 1].
      %   nmat_mesh           - (int array) The number of materials in each cell
      %                         of the mesh, size: [3, 3].
      %   matid_mesh          - (cell array of int arrays) The material IDs in each
      %                         cell of the mesh, with each cell containing an
      %                         int array of size [nmat_cell, 1].
      %   vf_mesh             - (cell array of double arrays) The volume fractions
      %                         of each material in each cell of the mesh, with each
      %                         cell containing a double array of size [nmat_cell, 1].
      %
      % Outputs:
      %   coords_final        - (double array) The final coordinates of the nodes
      %                         after reconstruction, size: [2*nnode_final, 1].
      %   nnode_final         - (int) The total number of nodes in the final
      %                         reconstructed geometry.
      %   nodes_for_interface - (int array) The nodes that lie on the material
      %                         interfaces, size: [nmat, 2].
      %   nodelist_for_mpoly  - (cell array of int arrays) The list of node indices
      %                         for each material polygon, each element of the cell
      %                         array is an int array.
      %   nnode_for_interface - (int array) The number of nodes for each material
      %                         interface, size: [nmat, 1].
      %   nnode_for_mpoly     - (int array) The number of nodes for each material
      %                         polygon, size: [nmat, 1].
      %
      % Original C declaration:
      % void reconstruct2d_nmat_pagosa(int geop,
      %                        double *xl, double *dx,
      %                        int **nmat_mesh, int ***matid_mesh, double ***vf_mesh,
      %                        int *nnode_final, double **coords_final,
      %                        int *nnode_for_interface, int **nodes_for_interface,
      %                        int *nnode_for_mpoly, int ***nodelist_for_mpoly)

      global util mesh;

      % Initialize variables
      dim = 2;
      xl = xl(1:dim);
      dx = dx(1:dim);
      dx_mx = max(dx);
      cell_vol = prod(dx ./ dx_mx);
      nnode = 4;
      coords = zeros(1, nnode*dim);

      % Define the initial coordinates for the 2D cell
      coords(1:2) = [  0.0,   0.0];
      coords(3:4) = [dx(1),   0.0] ./ dx_mx;
      coords(5:6) = [dx(1), dx(2)] ./ dx_mx;
      coords(7:8) = [  0.0, dx(2)] ./ dx_mx;

      nmat = nmat_mesh(2, 2);
      assert(nmat > 1, 'Number of materials must be greater than 1.');

      % Initialize output arrays
      coords_final = zeros((nnode + 2*nmat)*dim, 1);
      coords_final(1:nnode * dim) = coords(:);
      nnode_final = nnode;

      nnode_for_interface = zeros(nmat - 1, 1);
      nodes_for_interface = zeros(nmat - 1, 2);
      nnode_for_mpoly = zeros(nmat, 1);
      nodelist_for_mpoly = cell(nmat, 1);

      lsize = 0;

      % Initialize the volume fractions
      vfs = zeros(2, 3, 3);

      for m = 1:nmat-1
        % Calculate vfs for the current material
        vfs(1, 2, 2) = vfs(1, 2, 2) + vf_mesh(2, 2 ,m);
        vfs(2, 2, 2) = 1.0 - vfs(1, 2, 2);

        ids = matid_mesh(2, 2,:);
        for j = 1:3
          for i = 1:3
            if i == 2 && j == 2, continue, end
            nm_neighb = nmat_mesh(j, i);
            ids_neighb = matid_mesh(j, i,:);

            for m_neighb = 1:nm_neighb
              if ids_neighb(m_neighb) == ids(m)
                vfs(1, j, i) = vfs(1, j, i) + vf_mesh(j, i,m_neighb);
              else
                for m_other = m + 1:nmat
                  if ids_neighb(m_neighb) == ids(m_other)
                    vfs(2, j, i) = vfs(2, j, i) + vf_mesh(j, i,m_neighb);
                    break;
                  end
                end
              end
            end
          end
        end

        %--  now, we have the distribution of vfs[0] and vfs[1] each of which is a 3x3 mesh

        % Calculate the gradient and determine the priority
        priority_max = 0.0;
        idx = 1;
        norms = zeros(2, 2);

        for m_other = 1:2
          % mynorm => norms(m_other,:);
          grad = mesh.cal_cell_zgrad2d(dx, squeeze(vfs(m_other, :, :)));
          vsq = sum(grad .^ 2);
          factor = sqrt(vsq);
          priority = factor * sqrt(vfs(m_other, 2, 2));
          if priority > priority_max
            idx = m_other;
          end
          norms(m_other, :) = -grad / factor;
        end

        % Set the normal vector based on priority
        if idx == 2
          norm = -norms(idx, :);
        else
          norm = norms(idx, :);
        end

        % Call cal_distance2d method
        [~, nnode_new, coords_new, nnode_interface, nodelist_interface, ...
          nnode_lower, nodelist_lower, nnode_upper, nodelist_upper] = ...
          obj.cal_distance2d(geop, vfs(1, 2, 2), cell_vol, norm, nnode, coords, (1:4)');

        idx = nnode_final*dim;
        coords_final(idx+1:idx+nnode_new*dim) = coords_new(1:nnode_new*dim);

        nnode_for_interface(m) = nnode_interface;
        if nnode_interface == 2
          for i = 1:2
            n = nodelist_interface(i);
            if n <= nnode
              nodes_for_interface(m, i) = n;
            else
              nnew = n - nnode;
              nodes_for_interface(m, i) = nnode_final + nnew;
            end
          end
        end

        nnode_for_mpoly(m) = nnode_lower;
        for i = 1:nnode_lower
          n = nodelist_lower(i);
          if n <= nnode
            nodelist_for_mpoly{m}(i) = n;
          else
            nnew = n - nnode;
            nodelist_for_mpoly{m}(i) = nnode_final + nnew;
          end
        end

        if m < nmat - 1
          if ~isempty(nodelist_upper)
            nodelist_upper = [];
          end
        end
        nnode_previous = nnode_final;  % for the last nodelist_upper
        nnode_final = nnode_final + nnode_new;

      end % m-loop

      m = nmat;
      nnode_for_mpoly(m) = nnode_upper;
      for i = 1:nnode_upper
        n = nodelist_upper(i);
        if n <= nnode
          nodelist_for_mpoly{m}(i) = n;
        else
          nnew = n - nnode;
          nodelist_for_mpoly{m}(i) = nnode_previous + nnew;
        end
      end

      % Scale coords back
      for n = 1:nnode_final
        coords_final(1 + (n-1)*dim:n*dim) = xl' + dx_mx * coords_final(1 + (n-1)*dim:n*dim);
      end

      % Cleanup
      coords = [];
      coords_final = coords_final(1:nnode_final*2);
    end

    function [nnodex, coordsx] = remap2d_scaled(obj, ifdump, nnode, coords, inward_norm_ea_face, nnode_m, coords_m)
      % remap2d_scaled Scales, remaps, and rescales 2D coordinates for polygon intersection.
      %
      % This function scales the coordinates of two polygons, remaps their
      % intersection, and then rescales the coordinates back to their original
      % dimensions. The remapping is performed based on inward normals and
      % polygon geometry.
      %
      % Inputs:
      %   obj                - (c_vofm2d) The current object instance for VOF (Volume of Fluid) methods.
      %   ifdump             - (int) Indicator for whether to dump intermediate data (0 or 1).
      %   nnode              - (int) Number of nodes for the first polygon.
      %   coords             - (array of double) Coordinates of the first polygon [size: nnode x 2].
      %   inward_norm_ea_face - (array of double) Inward normals for the edges of the polygon [size: 4 x 2].
      %   nnode_m            - (int) Number of nodes for the second polygon.
      %   coords_m           - (array of double) Coordinates of the second polygon [size: nnode_m x 2].
      %
      % Outputs:
      %   nnodex             - (int) Number of nodes for the remapped intersecting polygon.
      %   coordsx            - (array of double) Coordinates of the remapped polygon [size: nnodex x 2].
      %
      % Uses:
      %   None.
      %
      % Modifies:
      %   None.
      %
      % Original C declaration:
      % void remap2d_scaled(int ifdump, int nnode, double *coords, double *inward_norm_ea_face,
      %                     int nnode_m, double *coords_m, int *nnodex, double **coordsx)

      dim = 2;

      % Allocate memory for scaled coordinates
      coords_scaled = zeros(nnode + nnode_m, dim);
      coords_m_scaled = coords_scaled(nnode+1:end, :);

      % Calculate maximum distance (dx_mx) for scaling
      c0 = coords(1, :);
      dx_mx = max(max(abs(coords(2:end, :) - c0), [], 2));

      % Shift and scale the first polygon
      coords_scaled(1:nnode, :) = (coords - c0) / dx_mx;

      % Shift and scale the second polygon
      coords_m_scaled(1:nnode_m, :) = (coords_m - c0) / dx_mx;

      % Call remap function to compute the intersection
      [nnodex, coordsx] = obj.remap2d_mat(ifdump, nnode, coords_scaled(1:nnode, :), ...
        inward_norm_ea_face, nnode_m, coords_m_scaled);

      % Rescale back to the original dimensions
      coordsx = c0 + dx_mx * coordsx;
    end

    function [nnodex, coordsx] = remap2d_mat(obj, ifdump, nnode, coords, inward_norm_ea_face, nnode_m, coords_m)
      % remap2d_mat Remaps 2D coordinates for polygon intersections.
      %
      % This function remaps two polygons based on their inward normals and
      % computes the intersection, returning the coordinates of the intersecting
      % polygon.
      %
      % Inputs:
      %   obj                - (c_vofm2d) The current object instance for VOF (Volume of Fluid) methods.
      %   ifdump             - (int) Indicator for whether to dump intermediate data (0 or 1).
      %   nnode              - (int) Number of nodes for the first polygon.
      %   coords             - (array of double) Coordinates of the first polygon [size: nnode x 2].
      %   inward_norm_ea_face - (array of double) Inward normals for the edges of the polygon [size: 4 x 2].
      %   nnode_m            - (int) Number of nodes for the second polygon.
      %   coords_m           - (array of double) Coordinates of the second polygon [size: nnode_m x 2].
      %
      % Outputs:
      %   nnodex             - (int) Number of nodes for the remapped intersecting polygon.
      %   coordsx            - (array of double) Coordinates of the remapped polygon [size: nnodex x 2].
      %
      % Original C declaration:
      % void remap2d_mat(int ifdump, int nnode, double *coords, double *inward_norm_ea_face,
      %                  int nnode_m, double *coords_m, int *nnodex, double **coordsx)

      dim = 2;
      accuracy = 1.0e-06;  % Accuracy threshold for being on the plane

      % Allocate memory for node lists and working arrays
      nnode_mx = max(nnode, nnode_m);
      n = 3 * nnode_mx;
      nodelist_default = 1:n;
      ds_ea_node = zeros(n, 1);
      coords_work = zeros(n, dim);
      coords_out = zeros(n, dim);
      coords_new = zeros(n, dim);

      my_nnode = nnode_m;
      my_coords = coords_m;
      nodelist_lower = [];
      nodelist_upper = [];

      % Iterate over the faces of the first polygon
      for n = 1:nnode
        n1 = mod(n, nnode) + 1;
        norm = inward_norm_ea_face(n, :);

        c = coords(n, :);
        c1 = coords(n1, :);
        distance = 0.5 * (norm(1) * (c(1) + c1(1)) + norm(2) * (c(2) + c1(2)));

        % Calculate the dot product for each node in the second polygon
        for i = 1:my_nnode
          ds_ea_node(i) = dot(norm, my_coords(i, :));
        end

        % Classify nodes based on their position relative to the plane
        node_loc = zeros(my_nnode, 1);
        node_loc(abs(ds_ea_node - distance) <= accuracy) = 0;  % On the plane
        node_loc(ds_ea_node > distance) = 1;  % Above the plane
        node_loc(ds_ea_node < distance) = -1;  % Below the plane

        % Find the interface and split the second polygon
        [nnode_new, coords_new, nnode_interface, nodes_interface, nnode_lower, nodelist_lower, nnode_upper, nodelist_upper] = ...
          obj.find_interface2d(my_nnode, my_coords, nodelist_default(1:my_nnode), norm, distance, node_loc);

        % Store the coordinates of the upper part of the intersecting polygon
        if nnode_upper < 3
          my_nnode = 0;
          break;
        else
          for i = 1:nnode_upper
            n0 = nodelist_upper(i);
            if n0 <= my_nnode
              coords_out(i, :) = my_coords(n0, :);
            else
              coords_out(i, :) = coords_new(n0 - my_nnode, :);
            end
          end
          my_nnode = nnode_upper;
          coords_work(1:my_nnode, :) = coords_out(1:my_nnode, :);
          my_coords = coords_work(1:my_nnode, :);
        end

        % Free temporary node lists
        nodelist_lower = [];
        nodelist_upper = [];
      end

      % Set output number of nodes and coordinates
      nnodex = my_nnode;
      if my_nnode >= 3
        coordsx = my_coords(1:my_nnode, :);
      else
        coordsx = [];
      end
    end

  end
end