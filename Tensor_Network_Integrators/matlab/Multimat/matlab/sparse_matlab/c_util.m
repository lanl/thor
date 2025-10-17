classdef c_util < handle
    % Utilities
    methods

%%%%%%%%%%%%%%% implemented and tested
        function area = cal_poly_area(~, nnode, coords, nnode_poly, nodelist)
            % cal_poly_area - Calculate the area of a polygon.
            %
            % This function calculates the area of a polygon defined by its vertices. The
            % polygon can be described directly by the coordinates of its nodes or by
            % referencing the nodes through an index list.
            %
            % Inputs:
            %   nnode       - (int) The total number of nodes available in the `coords` array.
            %   coords      - (double array) The coordinates of the polygon nodes,
            %                 size: [2*nnode, 1].
            %   nnode_poly  - (int) The number of nodes that define the polygon.
            %   nodelist    - (int array) A list of node indices that define the polygon,
            %                 size: [nnode_poly, 1]. If empty, the nodes are assumed to be
            %                 defined directly by their order in `coords`.
            %
            % Outputs:
            %   area        - (double) The computed area of the polygon.
            %
            % Original C declaration:
            % void cal_poly_area(int nnode, double *coords, int nnode_poly,
            %                    int *nodelist, double *area)
            % Initialize the area to 0
            area = 0.0;

            if ~isempty(nodelist)
                % Case when nodelist is provided
                for k = 1:nnode_poly
                    k1 = mod(k, nnode_poly) + 1;
                    n = nodelist(k);
                    n1 = nodelist(k1);
                    p0 = coords((n * 2 - 1):(n * 2));
                    p1 = coords((n1 * 2 - 1):(n1 * 2));
                    darea = p0(1) * p1(2) - p1(1) * p0(2);
                    area = area + darea;
                end
            else
                % Case when nodelist is not provided, using direct node order
                for k = 1:nnode_poly
                    k1 = mod(k, nnode_poly) + 1;
                    p0 = coords((k * 2 - 1):(k * 2));
                    p1 = coords((k1 * 2 - 1):(k1 * 2));
                    darea = p0(1) * p1(2) - p1(1) * p0(2);
                    area = area + darea;
                end
            end

            % Finalize the area calculation
            area = 0.5 * abs(area);
        end

        function vol = rz_area(~, nn, rz)
            % rz_area - Calculate the area under a curve in the rz-plane.
            %
            % This function calculates the volume (area under the curve) formed by a
            % polygon defined by points in the rz-plane. The points are assumed to
            % describe a closed polygon.
            %
            % Inputs:
            %   nn    - (int) The number of points defining the polygon in the rz-plane.
            %   rz    - (double array) The coordinates of the points in the rz-plane,
            %           given as a flattened array where each point has two values (r, z).
            %           Size: [2*nn, 1].
            %
            % Outputs:
            %   vol   - (double) The calculated volume (area under the curve) formed by the
            %           polygon.
            %
            % Original C declaration:
            % int rz_area(int nn, double *rz, double *vol)
            % Initialize volume to 0
            vol = 0.0;

            % Loop through the points to calculate the volume
            for k = 1:nn
                k1 = mod(k, nn) + 1;

                r0 = rz(2 * k - 1);
                z0 = rz(2 * k);
                r1 = rz(2 * k1 - 1);
                z1 = rz(2 * k1);

                % Calculate differential volume
                dz = z1 - z0;
                r02 = r0 * r0;
                r12 = r1 * r1;
                r0r1 = r0 * r1;

                dv = dz * (r12 + r0r1 + r02);
                vol = vol + dv;
            end

            % Finalize the volume calculation
            vol = vol / 3.0;
            vol = abs(vol);
        end

        function [node_order, ds_ea_node] = order_nodes_along_norm(~, dim, norm, ...
                                                                   nnode, coords)
            % order_nodes_along_norm - Order nodes along a specified normal vector.
            %
            % This function orders the nodes in a polygon based on their projection
            % along a given normal vector. The nodes are sorted according to their
            % distance along this vector, which is calculated as the dot product
            % between the node coordinates and the normal vector.
            %
            % Inputs:
            %   dim         - (int) The dimensionality of the space (e.g., 2 for 2D).
            %   norm        - (double array) The normal vector defining the direction
            %                 along which to order the nodes, size: [dim, 1].
            %   nnode       - (int) The number of nodes.
            %   coords      - (double array) The coordinates of the nodes, provided as a
            %                 flattened array where each point has `dim` values,
            %                 size: [dim*nnode, 1].
            %
            % Outputs:
            %   node_order  - (int array) The order of nodes based on their projection
            %                 along the normal vector, size: [nnode, 1].
            %   ds_ea_node  - (double array) The dot product of each node's coordinates
            %                 with the normal vector, size: [nnode, 1].
            %
            % Original C declaration:
            % void order_nodes_along_norm(int dim, double *norm,
            %                             int nnode, double *coords,
            %                             int *node_order, double *ds_ea_node)
            % Initialize outputs
            ds_ea_node = zeros(nnode, 1);
            node_order = zeros(nnode, 1);

            % Calculate the dot product of each node's coordinates with the normal vector
            for i = 1:nnode
                c = coords((i - 1) * dim + 1:i * dim);
                ds_ea_node(i) = sum(norm(1:dim) .* c);
            end

            % Create a list of pairs (ds_ea_node, index) and sort them
            [~, sorted_indices] = sort(ds_ea_node);

            % Update node_order based on the sorted indices
            node_order = sorted_indices;
        end
    
        function [sorted, sort_order] = r8sort(~, psi)
            % r8sort - Sort an array of double-precision numbers and track original indices.
            %
            % This function sorts an array of double-precision numbers (`psi`) in ascending 
            % order while preserving the original indices. The mapping of sorted indices 
            % back to the original order is stored in `order_tar_to_src`.
            %
            % Inputs:
            %   psi              - (double array) Array of numbers to be sorted, size: [n, 1].
            %
            % Outputs:
            %   sorted           - (double array) Sorted array of elements
            %   sort_order       - (int array) Mapping of sorted indices back to the original
            %                      order, size: [n, 1].
            % Original C declaration:
            % void r8sort(double *psi, int *order_tar_to_src, int n)
        
            [sorted, order_tar_to_src] = sort(psi, 'ascend');
            
        end

    end
end
