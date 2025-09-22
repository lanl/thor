classdef c_vof3d < handle
    
    properties
        accuracy = 1.0e-10;
    end
    methods

        function coords_scaled = coords_scale(~, ndim, x0_scale, dx_scale, nnode, coords)
            % coords_scale - Scale a set of coordinates by a specified factor.
            %
            % This function scales a set of `coords` based on a scaling factor (`dx_scale`) 
            % and an origin point (`x0_scale`). The resulting scaled coordinates are stored 
            % in `coords_scaled`.
            %
            % Inputs:
            %   ndim           - (int) Number of dimensions for each coordinate (e.g., 2 or 3).
            %   x0_scale       - (double array) Origin for scaling, size: [ndim, 1].
            %   dx_scale       - (double) Scaling factor.
            %   nnode          - (int) Number of nodes (coordinates) to scale.
            %   coords         - (double array) Original coordinates, size: [nnode, ndim].
            %
            % Outputs:
            %   coords_scaled  - (double array) Scaled coordinates, size: [nnode, ndim].
        
            % Calculate scaling factor
            factor = 1.0 / dx_scale;
        
            % Allocate output array for scaled coordinates
            coords_scaled = zeros(nnode, ndim);
        
            % Apply scaling to each coordinate
            for n = 1:nnode
                for i = 1:ndim
                    coords_scaled(n, i) = (coords(n, i) - x0_scale(i)) * factor;
                end
            end
        end

        function coords = coords_scale_back(~, ndim, x0_scale, dx_scale, nnode, coords)
            % coords_scale_back - Reverse the scaling of a set of coordinates.
            %
            % This function reverses a previous scaling transformation on a set of 
            % `coords` by applying the inverse of the scaling factor (`dx_scale`) 
            % and adding back the origin point (`x0_scale`). The coordinates are 
            % modified in place.
            %
            % Inputs:
            %   ndim        - (int) Number of dimensions for each coordinate (e.g., 2 or 3).
            %   x0_scale    - (double array) Origin for scaling, size: [ndim, 1].
            %   dx_scale    - (double) Scaling factor.
            %   nnode       - (int) Number of nodes (coordinates) to scale back.
            %   coords      - (double array) Coordinates to be scaled back, size: [nnode, ndim].
            %
            % Outputs:
            %   coords      - (double array) Coordinates after scaling back, modified in place.
        
            % Apply reverse scaling to each coordinate
            for n = 1:nnode
                for i = 1:ndim
                    coords(n, i) = dx_scale * coords(n, i) + x0_scale(i);
                end
            end
        end

        function coords_t = rotate_to_norm(~, norm, nnode, coords_s)
            % rotate_to_norm - Apply a rotation to align coordinates with a specified normal.
            %
            % This function rotates a set of 3D coordinates (`coords_s`) so that the z-axis
            % aligns with the specified normal vector (`norm`) and returns the transformed 
            % coordinates (`coords_t`).
            %
            % Inputs:
            %   norm       - (double array) Normal vector defining the target z-axis orientation,
            %                size: [3, 1].
            %   nnode      - (int) Number of nodes defining the polygon.
            %   coords_s   - (double array) Source coordinates of nodes in the original frame,
            %                size: [nnode, 3].
            %
            % Outputs:
            %   coords_t   - (double array) Transformed coordinates after aligning with `norm`,
            %                size: [nnode, 3].
        
            global small;
        
            % Initialize rotation parameters
            if abs(norm(3) - 1.0) < small  % No change needed
                cost = 1.0; sint = 0.0; sinp = 0.0; cosp = 1.0;
            elseif abs(norm(3) + 1.0) < small  % -z -> z, y -> x, x -> y
                cost = -1.0; sint = 0.0; sinp = -1.0; cosp = 0.0;
            elseif abs(norm(1) - 1.0) < small  % x -> z, -z -> x, y -> y
                cost = 0.0; sint = 1.0; sinp = 0.0; cosp = 1.0;
            elseif abs(norm(1) + 1.0) < small  % -x -> z, -y -> y, -z -> x
                cost = 0.0; sint = 1.0; sinp = 0.0; cosp = -1.0;
            elseif abs(norm(2) - 1.0) < small  % y -> z, -z -> x, -x -> y
                cost = 0.0; sint = 1.0; sinp = 1.0; cosp = 0.0;
            elseif abs(norm(2) + 1.0) < small  % -y -> z, -z -> x, x -> y
                cost = 0.0; sint = 1.0; sinp = -1.0; cosp = 0.0;
            else  % General case
                cost = norm(3);
                sint = sqrt(norm(1)^2 + norm(2)^2);
                if sint > 0.0
                    tmp = 1.0 / sint;
                    sinp = norm(2) * tmp;
                    cosp = norm(1) * tmp;
                else
                    sinp = 1.0;
                    cosp = 0.0;
                end
            end
        
            % Allocate output array for transformed coordinates
            coords_t = zeros(nnode, 3);
        
            % Apply rotation to each node
            for i = 1:nnode
                c = coords_s(i, :);  % Current coordinates
        
                % Rotate coordinates to align with the normal vector
                coords_t(i, 1) =  cost * cosp * c(1) + cost * sinp * c(2) - sint * c(3);
                coords_t(i, 2) = -      sinp * c(1) +      cosp * c(2);
                coords_t(i, 3) =  sint * cosp * c(1) + sint * sinp * c(2) + cost * c(3);
            end
        end

        function coords_t = rotate_back(~, norm, nnode, coords_s)
            % rotate_back - Apply a rotation to transform coordinates back to a standard frame.
            %
            % This function rotates a set of 3D coordinates (`coords_s`) based on a 
            % normal vector (`norm`) and returns the transformed coordinates (`coords_t`).
            % It adjusts the rotation angles depending on the orientation of `norm`.
            %
            % Inputs:
            %   norm       - (double array) Normal vector defining the orientation 
            %                of the interface, size: [3, 1].
            %   nnode      - (int) Number of nodes defining the polygon.
            %   coords_s   - (double array) Source coordinates of nodes in the current 
            %                frame, size: [nnode, 3].
            %
            % Outputs:
            %   coords_t   - (double array) Transformed coordinates after applying 
            %                the rotation, size: [nnode, 3].
            
            global small;
            
            % Initialize rotation parameters
            if abs(norm(3) - 1.0) < small  % no change
                cost = 1.0; sint = 0.0; sinp = 0.0; cosp = 1.0;
            elseif abs(norm(3) + 1.0) < small  % -z -> z, y -> x, x -> y
                cost = -1.0; sint = 0.0; sinp = -1.0; cosp = 0.0;
            elseif abs(norm(1) - 1.0) < small  % x -> z, -z -> x, y -> y
                cost = 0.0; sint = 1.0; sinp = 0.0; cosp = 1.0;
            elseif abs(norm(1) + 1.0) < small  % -x -> z, -y -> y, -z -> x
                cost = 0.0; sint = 1.0; sinp = 0.0; cosp = -1.0;
            elseif abs(norm(2) - 1.0) < small  % y -> z, -z -> x, -x -> y
                cost = 0.0; sint = 1.0; sinp = 1.0; cosp = 0.0;
            elseif abs(norm(2) + 1.0) < small  % -y -> z, -z -> x, x -> y
                cost = 0.0; sint = 1.0; sinp = -1.0; cosp = 0.0;
            else  % General case
                cost = norm(3);
                sint = sqrt(norm(1)^2 + norm(2)^2);
                if sint > 0.0
                    tmp = 1.0 / sint;
                    sinp = norm(2) * tmp;
                    cosp = norm(1) * tmp;
                else
                    sinp = 1.0;
                    cosp = 0.0;
                end
            end
        
            % Allocate output array for transformed coordinates
            coords_t = zeros(nnode, 3);
        
            % Apply rotation to each node
            for i = 1:nnode
                c = coords_s(i, :);  % Current coordinates
        
                % Rotate coordinates
                coords_t(i, 1) =  cost*cosp*c(1) - sinp*c(2) + sint*cosp*c(3);
                coords_t(i, 2) =  cost*sinp*c(1) + cosp*c(2) + sint*sinp*c(3);
                coords_t(i, 3) = -sint*c(1)                       + cost*c(3);
            end
        end

        function [nedge, edgelist_for_face, nodelist_for_edge] = get_edges(~, ...
                  nface, nnode, nnode_for_face, nodelist_for_face)
            % get_edges - Determine edges for each face and avoid duplicates.
            %
            % This function generates a unique edge list and identifies edges for 
            % each face in a mesh, avoiding duplicate edges shared by adjacent faces.
            %
            % Inputs:
            %   nface               - (int) Number of faces in the mesh.
            %   nnode               - (int) Number of nodes in the mesh.
            %   nnode_for_face      - (int array) Number of nodes per face, size: [nface, 1].
            %   nodelist_for_face   - (int array) List of node indices for each face, size: [sum(nnode_for_face), 1].
            %
            % Outputs:
            %   nedge               - (int) Number of unique edges in the mesh.
            %   edgelist_for_face   - (int array) Edge indices for each face, size: [sum(nnode_for_face), 1].
            %   nodelist_for_edge   - (int array) Node pairs for each unique edge, size: [nedge, 2].
        
            % Initialize variables
            ne = sum(nnode_for_face);  % Total number of edges, assuming max (one per node)
            nodelist_for_edge_work = zeros(2 * ne, 1);  % Temporary storage for edges
            edgelist_for_face = zeros(ne, 1);           % Edge indices for each face
        
            % Pointers for easy manipulation of nodelist and edge list
            nodelist_s = nodelist_for_face;
            nodelist_t = nodelist_for_edge_work;
            edgelist_t = edgelist_for_face;
        
            nedge = 0;  % Initialize unique edge count
            edge_pos = 1;  % Position in `edgelist_for_face`
        
            % Iterate over each face
            for f = 1:nface
                nn = nnode_for_face(f);  % Number of nodes in the current face
        
                % Iterate over each edge in the face
                for i = 1:nn
                    i1 = mod(i, nn) + 1;  % Wrap-around to first node
                    n0 = nodelist_s(i);
                    n1 = nodelist_s(i1);
        
                    % Check if this edge already exists
                    e = -1;
                    for k = 1:nedge
                        m0 = nodelist_for_edge_work(2 * (k - 1) + 1);
                        m1 = nodelist_for_edge_work(2 * (k - 1) + 2);
                        if ((m0 == n0 && m1 == n1) || (m0 == n1 && m1 == n0))
                            e = k - 1;  % Adjust for 0-based indexing
                            break;
                        end
                    end
        
                    % Assign edge index, either existing or new
                    if e >= 0
                        edgelist_t(edge_pos) = e;
                    else
                        edgelist_t(edge_pos) = nedge;
                        nodelist_t(2 * nedge + 1) = n0;
                        nodelist_t(2 * nedge + 2) = n1;
                        nedge = nedge + 1;
                    end
                    edge_pos = edge_pos + 1;
                end
                nodelist_s = nodelist_s(nn + 1:end);  % Move to next face nodes
            end
        
            % Finalize the unique edge list
            nodelist_for_edge = reshape(nodelist_for_edge_work(1:2 * nedge), [nedge, 2]);
        end

        function [facelist_for_edge, nedge_ea_node, edgelist_for_node] = get_faces_ea_edge(~, ...
                nface, nedge, nnode, nnode_ea_face, nodelist_for_face, edgelist_for_face, nodelist_for_edge)
            % get_faces_ea_edge - Determine faces associated with each edge and edges associated with each node.
            %
            % This function determines which faces are associated with each edge and 
            % which edges are associated with each node in a 3D mesh. It creates a list 
            % of edges for each node and a list of faces for each edge.
            %
            % Inputs:
            %   nface              - (int) Number of faces in the mesh.
            %   nedge              - (int) Number of unique edges in the mesh.
            %   nnode              - (int) Number of nodes in the mesh.
            %   nnode_ea_face      - (int array) Number of nodes per face, size: [nface, 1].
            %   nodelist_for_face  - (int array) Node indices for each face, size: [sum(nnode_ea_face), 1].
            %   edgelist_for_face  - (int array) Edge indices for each face, size: [sum(nnode_ea_face), 1].
            %   nodelist_for_edge  - (int array) Node pairs for each unique edge, size: [nedge, 2].
            %
            % Outputs:
            %   facelist_for_edge  - (int array) Face indices for each edge (2 per edge), size: [2*nedge, 1].
            %   nedge_ea_node      - (int array) Number of edges associated with each node, size: [nnode, 1].
            %   edgelist_for_node  - (int array) List of edges associated with each node, size: [sum(nedge_ea_node), 1].
        
            % Initialize facelist_for_edge with -1
            facelist_for_edge = -1 * ones(2 * nedge, 1);
        
            % Populate facelist_for_edge by associating faces with edges
            offset = 1;
            for f = 1:nface
                ne = nnode_ea_face(f);  % Number of edges for this face
                edges = edgelist_for_face(offset:offset + ne - 1);  % Edges for this face
                for i = 1:ne
                    e = edges(i);        % Current edge index
                    e2 = 2 * (e - 1) + 1;  % Index in facelist_for_edge (MATLAB 1-based indexing)
        
                    % Assign face index to the appropriate position in facelist_for_edge
                    if facelist_for_edge(e2) == -1
                        facelist_for_edge(e2) = f;
                    else
                        facelist_for_edge(e2 + 1) = f;
                    end
                end
                offset = offset + ne;
            end
        
            % Initialize nedge_ea_node with zeros and calculate the size of edgelist_for_node
            nedge_ea_node = zeros(nnode, 1);
            sz = 0;
            for e = 1:nedge
                nodes = nodelist_for_edge(2 * (e - 1) + 1:2 * e);  % Node indices for the edge
                for i = 1:2
                    n = nodes(i);
                    nedge_ea_node(n) = nedge_ea_node(n) + 1;  % Increment edge count for this node
                    sz = sz + 1;
                end
            end
        
            % Allocate edgelist_for_node with the calculated size
            edgelist_for_node = zeros(sz, 1);
        
            % Create offset for each node in edgelist_for_node
            offset_ea_node = cumsum([0; nedge_ea_node(1:end - 1)]) + 1;  % Start positions for each node
        
            % Reset nedge_ea_node to track positions in edgelist_for_node
            nedge_ea_node(:) = 0;
        
            % Populate edgelist_for_node by associating edges with nodes
            for e = 1:nedge
                nodes = nodelist_for_edge(2 * (e - 1) + 1:2 * e);  % Node indices for the edge
                for i = 1:2
                    n = nodes(i);
                    pos = offset_ea_node(n) + nedge_ea_node(n);
                    edgelist_for_node(pos) = e;  % Assign edge to the node's edge list
                    nedge_ea_node(n) = nedge_ea_node(n) + 1;  % Increment edge count for this node
                end
            end
        end

        function interior_point = get_interior_point(~, dim, nface, nnode, ...
                 coords, nnode_ea_face, nodelist_for_face)
            % get_interior_point - Calculate an interior point for a 3D shape.
            %
            % This function calculates an interior point by averaging the coordinates 
            % of unique nodes that belong to any face in the mesh. The resulting 
            % interior point is the centroid of all the nodes associated with faces.
            %
            % Inputs:
            %   dim               - (int) Dimension of each coordinate (e.g., 3 for 3D).
            %   nface             - (int) Number of faces in the mesh.
            %   nnode             - (int) Number of nodes in the mesh.
            %   coords            - (double array) Coordinates of all nodes, size: [nnode, dim].
            %   nnode_ea_face     - (int array) Number of nodes for each face, size: [nface, 1].
            %   nodelist_for_face - (int array) Node indices for each face, size: [sum(nnode_ea_face), 1].
            %
            % Outputs:
            %   interior_point    - (double array) Calculated interior point, size: [1, dim].
        
            % Initialize `to_include` to mark nodes that belong to any face
            to_include = zeros(nnode, 1);
            offset = 1;
        
            % Mark nodes that belong to each face
            for f = 1:nface
                nn = nnode_ea_face(f);
                nodes = nodelist_for_face(offset:offset + nn - 1);
                to_include(nodes) = 1;
                offset = offset + nn;
            end
        
            % Initialize the interior point as a zero vector
            interior_point = zeros(1, dim);
        
            % Sum the coordinates of nodes marked in `to_include`
            nn = 0;  % Counter for the included nodes
            for i = 1:nnode
                if to_include(i) == 0
                    continue;
                end
                interior_point = interior_point + coords(i, :);
                nn = nn + 1;
            end
        
            % Average the coordinates to find the centroid (interior point)
            if nn > 0
                interior_point = interior_point / nn;
            else
                warning('No nodes found in the faces.');
            end
        end

        function [nodelist_new, node2edge_new, err] = sort_nodes(~, nnode, coords, node2edge)
            % sort_nodes - Sort nodes in a 2D plane based on their angle relative to a center point.
            %
            % This function calculates the center of a set of 2D nodes, then sorts the 
            % nodes based on their angle from the center relative to the x-axis. The 
            % sorted node indices and corresponding edges are returned.
            %
            % Inputs:
            %   nnode         - (int) Number of nodes.
            %   coords        - (double array) Coordinates of nodes, size: [nnode, 2].
            %   node2edge     - (int array) Edge indices associated with each node, size: [nnode, 1].
            %
            % Outputs:
            %   nodelist_new  - (int array) Sorted indices of nodes based on angle, size: [nnode, 1].
            %   node2edge_new - (int array) Edge indices sorted to match the sorted nodes, size: [nnode, 1].
            %   err           - (int) Error flag, 0 if no error occurs.
        
            global util;

            % Initialize outputs
            err = 0;
            pi2 = 2 * pi;
            dim = 2;
            
            % Step 1: Calculate the center point (ctr)
            ctr = mean(coords, 1);  % Average of the coordinates in x and y
        
            % Step 2: Calculate the angle from the center to each node
            a_ea_node = zeros(nnode, 1);
            for n = 1:nnode
                dx = ctr - coords(n, :);  % Vector from node to center
                cosa = dx(1) / norm(dx);  % Cosine of the angle
                a = acos(cosa);            % Angle in radians
                if dx(2) < 0
                    a = pi2 - a;
                end
                a_ea_node(n) = a;
            end
        
            % Step 3: Sort the nodes by their angle
            [~, tar_to_src] = sort(a_ea_node, 'ascend');  % Sort angles, get sorted indices
        
            % Step 4: Reorder nodelist and node2edge based on the sorted order
            nodelist_new = tar_to_src;           % Sorted indices of nodes
            node2edge_new = node2edge(tar_to_src);  % Edge indices sorted accordingly
        end

        function point = line_intersect_z(~, c0, c1, zint)
            % line_intersect_z - Find the intersection of a line segment with a given z-plane.
            %
            % This function calculates the intersection point between a line segment, 
            % defined by points `c0` and `c1`, and a plane with a constant z-value (`zint`).
            %
            % Inputs:
            %   c0     - (double array) Coordinates of the first point on the line, size: [1, 3].
            %   c1     - (double array) Coordinates of the second point on the line, size: [1, 3].
            %   zint   - (double) The z-value of the intersection plane.
            %
            % Outputs:
            %   point  - (double array) The intersection point on the z-plane, size: [1, 3].
        
            point = zeros(1, 3);
      
            % Interpolate the x and y coordinates of the intersection point
            t = (zint - c0(3)) / (c1(3) - c0(3));
            point(1:2) = c0(1:2) + (c1(1:2) - c0(1:2)) * t;
            point(3) = zint;
        end

        function [nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower, ...
                  nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper, ...
                  faceindexlow_of_upper_interface, faceindexhgh_of_upper_interface] = ...
                  get_lower_upper_polyhedron(~, nface, nedge, nnode, nnode_ea_face, nodelist_for_face, ...
                                             edgelist_for_face, nodelist_for_edge, nnode_new, edge_to_node, ...
                                             nnode_interface, nodes_interface, node_loc, edge_loc, ...
                                             interface_is_face)
            % get_lower_upper_polyhedron - Generate lower and upper polyhedrons from a 3D mesh.
            %
            % This function partitions a 3D mesh into lower and upper polyhedrons based on
            % specified interface nodes, edges, and faces. It classifies nodes, edges, and 
            % faces as belonging to either the lower or upper polyhedron, organizing them 
            % accordingly.
            %
            % Inputs:
            %   nface                             - (int) Total number of faces in the original mesh.
            %   nedge                             - (int) Total number of edges in the original mesh.
            %   nnode                             - (int) Total number of nodes in the original mesh.
            %   nnode_ea_face                     - (int array) Number of nodes per face, size: [nface, 1].
            %   nodelist_for_face                 - (int array) Node indices for each face, size: [sum(nnode_ea_face), 1].
            %   edgelist_for_face                 - (int array) Edge indices for each face, size: [sum(nnode_ea_face), 1].
            %   nodelist_for_edge                 - (int array) Node pairs for each edge, size: [nedge, 2].
            %   nnode_new                         - (int) Number of new nodes created on the interface.
            %   edge_to_node                      - (int array) Mapping from edges to nodes, size: [nedge, 2].
            %   nnode_interface                   - (int) Number of nodes on the interface.
            %   nodes_interface                   - (int array) Indices of nodes on the interface, size: [nnode_interface, 1].
            %   node_loc                          - (int array) Classification of nodes (-1 for lower, 
            %                                       0 for on interface, 1 for upper), size: [nnode, 1].
            %   edge_loc                          - (int array) Classification of edges (-1 for lower, 
            %                                       0 for on interface, 1 for upper), size: [nedge, 1].
            %   interface_is_face                 - (int) Flag indicating if the interface is treated 
            %                                       as a face (1 if true, 0 otherwise).
            %   faceindexlow_of_upper_interface   - (int) Index of the lower polyhedron face on the 
            %                                       upper interface.
            %   faceindexhgh_of_upper_interface   - (int) Index of the upper polyhedron face on the 
            %                                       upper interface.
            %
            % Outputs:
            %   nface_lower                       - (int) Number of faces in the lower polyhedron.
            %   nnode_ea_face_lower               - (int array) Number of nodes per face for the lower 
            %                                       polyhedron, size: [nface_lower, 1].
            %   nodelist_for_face_lower           - (int array) Node indices for each face in the lower 
            %                                       polyhedron, size: [sum(nnode_ea_face_lower), 1].
            %   face2old_lower                    - (int array) Mapping from lower polyhedron faces 
            %                                       to the original mesh faces, size: [nface_lower, 1].
            %   nface_upper                       - (int) Number of faces in the upper polyhedron.
            %   nnode_ea_face_upper               - (int array) Number of nodes per face for the upper 
            %                                       polyhedron, size: [nface_upper, 1].
            %   nodelist_for_face_upper           - (int array) Node indices for each face in the upper 
            %                                       polyhedron, size: [sum(nnode_ea_face_upper), 1].
            %   face2old_upper                    - (int array) Mapping from upper polyhedron faces 
            %                                       to the original mesh faces, size: [nface_upper, 1].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower - 
            %       Outputs defining the lower polyhedron.
            %   nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper - 
            %       Outputs defining the upper polyhedron.
            %
            % Original C declaration:
            % void get_lower_upper_polyhedron(int nface, int nedge, int nnode,
            %                                 int *nnode_ea_face, int *nodelist_for_face,
            %                                 int *edgelist_for_face, int *nodelist_for_edge,
            %                                 int nnode_new, int *edge_to_node,
            %                                 int nnode_interface, int *nodes_interface,
            %                                 int *node_loc, int *edge_loc,
            %                                 int interface_is_face,
            %              int *faceindexlow_of_upper_interface, int *faceindexhgh_of_upper_interface,
            %              int *nface_lower, int *nnode_ea_face_lower, int *nodelist_for_face_lower, int *face2old_lower,
            %              int *nface_upper, int *nnode_ea_face_upper, int *nodelist_for_face_upper, int *face2old_upper);
        
            % Initialization
            faceindexlow_of_upper_interface = -1;
            faceindexhgh_of_upper_interface = -1;
            edge_to_include = zeros(nedge, 1);
            face_to_include = zeros(nface, 1);
            face_all_in = zeros(nface, 1);
        
            % Determine edges to include in lower polyhedron
            for e = 1:nedge
                nodes = nodelist_for_edge(2 * (e - 1) + (1:2));
                n0 = nodes(1);
                n1 = nodes(2);
                loc0 = node_loc(n0);
                loc1 = node_loc(n1);
        
                if ((loc0 <= 0 && loc1 <= 0) || (loc0 < 0 && loc1 >= 0) || (loc0 >= 0 && loc1 < 0))
                    edge_to_include(e) = 1;
                end
            end
        
            % Determine faces to include in lower polyhedron
            offset = 1;
            for f = 1:nface
                nn = nnode_ea_face(f);
                nodes = nodelist_for_face(offset:offset + nn - 1);
                some_node_inside = any(node_loc(nodes) < 0);
                some_node_outside = any(node_loc(nodes) > 0);
        
                if some_node_inside && ~some_node_outside
                    face_all_in(f) = 1;
                    face_to_include(f) = 1;
                elseif ~some_node_inside && some_node_outside
                    face_to_include(f) = 0;
                    face_all_in(f) = 0;
                elseif some_node_inside && some_node_outside
                    face_to_include(f) = 1;
                    face_all_in(f) = 0;
                end
                offset = offset + nn;
            end
        
            % Initialize lower polyhedron
            nface_lower = 0;
            offset_t = 1;
        
            % Include interface in lower polyhedron if needed
            if ~interface_is_face && nnode_interface >= 3
                nnode_ea_face_lower(nface_lower + 1) = nnode_interface;
                nodelist_for_face_lower(offset_t:offset_t + nnode_interface - 1) = nodes_interface;
                faceindexlow_of_upper_interface = nface_lower;
                face2old_lower(nface_lower + 1) = -1;
                offset_t = offset_t + nnode_interface;
                nface_lower = nface_lower + 1;
            end
        
            % Process each face in the lower polyhedron
            offset = 1;
            for f = 1:nface
                ne = nnode_ea_face(f);
                if face_all_in(f)
                    nnode_ea_face_lower(nface_lower + 1) = ne;
                    nodelist_for_face_lower(offset_t:offset_t + ne - 1) = nodelist_for_face(offset:offset + ne - 1);
                    face2old_lower(nface_lower + 1) = f;
                    offset_t = offset_t + ne;
                    nface_lower = nface_lower + 1;
                elseif face_to_include(f)
                    nn = 0;
                    nodes_t = zeros(ne + nnode_interface, 1);
                    edges = edgelist_for_face(offset:offset + ne - 1);
                    nodes = nodelist_for_face(offset:offset + ne - 1);
        
                    % Process edges of each face
                    for i = 1:ne
                        i1 = mod(i, ne) + 1;
                        e = edges(i);
                        e1 = edges(i1);
                        nodes_on_e = nodelist_for_edge(2 * (e - 1) + (1:2));
                        nodes_on_e1 = nodelist_for_edge(2 * (e1 - 1) + (1:2));
                        n0 = nodes_on_e(1);
                        n1 = nodes_on_e(2);
        
                        if i == 1
                            if any(n0 == nodes_on_e1)
                                last_n = n1;
                            elseif any(n1 == nodes_on_e1)
                                last_n = n0;
                            else
                                error('ERROR: last_n');
                            end
                        end
        
                        if edge_to_node(e) >= 0
                            if node_loc(last_n) > 0
                                nodes_t(nn + 1) = edge_to_node(e);
                                nn = nn + 1;
                                nodes_t(nn + (1:2)) = nodes_on_e(node_loc(nodes_on_e) == -1);
                                nn = nn + sum(node_loc(nodes_on_e) == -1);
                            else
                                nodes_t(nn + (1:2)) = nodes_on_e(node_loc(nodes_on_e) == -1);
                                nn = nn + sum(node_loc(nodes_on_e) == -1);
                                nodes_t(nn + 1) = edge_to_node(e);
                                nn = nn + 1;
                            end
                        elseif edge_to_include(e)
                            if n0 == last_n
                                if node_loc(n0) <= 0, nodes_t(nn + 1) = n0; nn = nn + 1; end
                                if node_loc(n1) <= 0, nodes_t(nn + 1) = n1; nn = nn + 1; end
                            elseif n1 == last_n
                                if node_loc(n1) <= 0, nodes_t(nn + 1) = n1; nn = nn + 1; end
                                if node_loc(n0) <= 0, nodes_t(nn + 1) = n0; nn = nn + 1; end
                            end
                        end
                        last_n = nodes_on_e1(any(nodes_on_e1 == [n0, n1]));
                    end
        
                    % Consolidate and remove redundant nodes
                    if nn > 1 && nodes_t(nn) == nodes_t(1), nn = nn - 1; end
                    nodes_t = unique(nodes_t(1:nn), 'stable');
                    nn = length(nodes_t);
        
                    if nn >= 3
                        nnode_ea_face_lower(nface_lower + 1) = nn;
                        nodelist_for_face_lower(offset_t:offset_t + nn - 1) = nodes_t;
                        face2old_lower(nface_lower + 1) = f;
                        offset_t = offset_t + nn;
                        nface_lower = nface_lower + 1;
                    end
                end
                offset = offset + ne;
            end
        
            % Process upper polyhedron
            offset_t = 1;
            if ~interface_is_face && nnode_interface >= 3
                nnode_ea_face_upper(nface_upper + 1) = nnode_interface;
                nodelist_for_face_upper(offset_t:offset_t + nnode_interface - 1) = nodes_interface;
                faceindexhgh_of_upper_interface = nface_upper;
                face2old_upper(nface_upper + 1) = -1;
                offset_t = offset_t + nnode_interface;
                nface_upper = nface_upper + 1;
            end
        
            offset = 1;
            for f = 1:nface
                ne = nnode_ea_face(f);
                if face_all_in(f)
                    nnode_ea_face_upper(nface_upper + 1) = ne;
                    nodelist_for_face_upper(offset_t:offset_t + ne - 1) = nodelist_for_face(offset:offset + ne - 1);
                    face2old_upper(nface_upper + 1) = f;
                    offset_t = offset_t + ne;
                    nface_upper = nface_upper + 1;
                elseif face_to_include(f)
                    nn = 0;
                    nodes_t = zeros(ne + nnode_interface, 1);
                    edges = edgelist_for_face(offset:offset + ne - 1);
                    nodes = nodelist_for_face(offset:offset + ne - 1);
        
                    for i = 1:ne
                        i1 = mod(i, ne) + 1;
                        e = edges(i);
                        e1 = edges(i1);
                        nodes_on_e = nodelist_for_edge(2 * (e - 1) + (1:2));
                        nodes_on_e1 = nodelist_for_edge(2 * (e1 - 1) + (1:2));
                        n0 = nodes_on_e(1);
                        n1 = nodes_on_e(2);
        
                        if i == 1
                            if any(n0 == nodes_on_e1)
                                last_n = n1;
                            elseif any(n1 == nodes_on_e1)
                                last_n = n0;
                            else
                                error('ERROR: last_n');
                            end
                        end
        
                        if edge_to_node(e) >= 0
                            if node_loc(last_n) < 0
                                nodes_t(nn + 1) = edge_to_node(e);
                                nn = nn + 1;
                                nodes_t(nn + (1:2)) = nodes_on_e(node_loc(nodes_on_e) == 1);
                                nn = nn + sum(node_loc(nodes_on_e) == 1);
                            else
                                nodes_t(nn + (1:2)) = nodes_on_e(node_loc(nodes_on_e) == 1);
                                nn = nn + sum(node_loc(nodes_on_e) == 1);
                                nodes_t(nn + 1) = edge_to_node(e);
                                nn = nn + 1;
                            end
                        elseif edge_to_include(e)
                            if n0 == last_n
                                if node_loc(n0) >= 0, nodes_t(nn + 1) = n0; nn = nn + 1; end
                                if node_loc(n1) >= 0, nodes_t(nn + 1) = n1; nn = nn + 1; end
                            elseif n1 == last_n
                                if node_loc(n1) >= 0, nodes_t(nn + 1) = n1; nn = nn + 1; end
                                if node_loc(n0) >= 0, nodes_t(nn + 1) = n0; nn = nn + 1; end
                            end
                        end
                        last_n = nodes_on_e1(any(nodes_on_e1 == [n0, n1]));
                    end
        
                    if nn > 1 && nodes_t(nn) == nodes_t(1), nn = nn - 1; end
                    nodes_t = unique(nodes_t(1:nn), 'stable');
                    nn = length(nodes_t);
        
                    if nn >= 3
                        nnode_ea_face_upper(nface_upper + 1) = nn;
                        nodelist_for_face_upper(offset_t:offset_t + nn - 1) = nodes_t;
                        face2old_upper(nface_upper + 1) = f;
                        offset_t = offset_t + nn;
                        nface_upper = nface_upper + 1;
                    end
                end
                offset = offset + ne;
            end
        
            if interface_is_face
                if nface_lower > 3
                    faceindexlow_of_upper_interface = interface_is_face - 1;
                end
                if nface_upper > 3
                    faceindexhgh_of_upper_interface = interface_is_face - 1;
                end
            end
        
            % Free memory
            clear edge_to_include;
        end

        function vol = tet_vol(~, a, b, c, d)
            % tet_vol - Calculate the volume of a tetrahedron defined by four points.
            %
            % This function calculates the volume of a tetrahedron defined by points 
            % `a`, `b`, `c`, and `d` in 3D space using vector cross products.
            %
            % Inputs:
            %   a  - (double array) Coordinates of the first vertex, size: [1, 3].
            %   b  - (double array) Coordinates of the second vertex, size: [1, 3].
            %   c  - (double array) Coordinates of the third vertex, size: [1, 3].
            %   d  - (double array) Coordinates of the fourth vertex, size: [1, 3].
            %
            % Outputs:
            %   vol - (double) Calculated volume of the tetrahedron.
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   None.
            %
            % Original C declaration:
            % double tet_vol(double *a, double *b, double *c, double *d)
        
            % Compute vectors ad, bd, and cd relative to point d
            ad = a - d;
            bd = b - d;
            cd = c - d;
        
            % Calculate the cross product of bd and cd
            cross = [bd(2)*cd(3) - bd(3)*cd(2), ...
                     bd(3)*cd(1) - bd(1)*cd(3), ...
                     bd(1)*cd(2) - bd(2)*cd(1)];
        
            % Calculate the volume
            vol = abs(dot(ad, cross)) / 6.0;
        end

        function ctr = cal_zone_ctr(~, nface, nnode, coords, nnode_for_face, nodelist_for_face)
            % cal_zone_ctr - Calculate the centroid of a 3D zone defined by its faces.
            %
            % This function calculates the centroid of a 3D zone (e.g., polyhedron) 
            % by averaging the coordinates of the unique nodes that are part of the 
            % zone’s faces.
            %
            % Inputs:
            %   nface             - (int) Number of faces defining the zone.
            %   nnode             - (int) Total number of nodes in the mesh.
            %   coords            - (double array) Coordinates of all nodes, size: [nnode, 3].
            %   nnode_for_face    - (int array) Number of nodes for each face, size: [nface, 1].
            %   nodelist_for_face - (int array) Node indices for each face, size: [sum(nnode_for_face), 1].
            %
            % Outputs:
            %   ctr               - (double array) Calculated centroid of the zone, size: [1, 3].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   None.
            %
            % Original C declaration:
            % void cal_zone_ctr(int nface, int nnode, double *coords, int *nnode_for_face,
            %                   int *nodelist_for_face, double *ctr)
        
            % Initialize included nodes array to mark nodes in faces
            included = zeros(nnode, 1);
        
            % Mark nodes that are part of any face in the zone
            offset = 1;
            for f = 1:nface
                nn = nnode_for_face(f);
                nodes = nodelist_for_face(offset:offset + nn - 1);
                included(nodes) = 1;
                offset = offset + nn;
            end
        
            % Initialize centroid
            ctr = zeros(1, 3);
        
            % Sum coordinates of included nodes to calculate the centroid
            nn = 0;  % Counter for included nodes
            for n = 1:nnode
                if included(n)
                    ctr = ctr + coords(n, :);
                    nn = nn + 1;
                end
            end
        
            % Divide by the number of included nodes to get the average (centroid)
            if nn > 0
                ctr = ctr / nn;
            else
                warning('No nodes included in the zone.');
            end
        end

        function vol = cal_vol(obj, nface, nnode, coords, nnode_for_face, nodelist_for_face)
            % cal_vol - Calculate the volume of a polyhedral zone in 3D space.
            %
            % This function calculates the volume of a polyhedral zone defined by 
            % nodes and faces. For tetrahedrons, it calculates the volume directly; 
            % for other polyhedrons, it uses the centroid and face centroids.
            %
            % Inputs:
            %   obj               - (c_vof3d) Encapsulating class object.
            %   nface             - (int) Number of faces in the polyhedral zone.
            %   nnode             - (int) Number of nodes in the polyhedral zone.
            %   coords            - (double array) Coordinates of all nodes, size: [nnode, 3].
            %   nnode_for_face    - (int array) Number of nodes per face, size: [nface, 1].
            %   nodelist_for_face - (int array) Node indices for each face, 
            %                       size: [sum(nnode_for_face), 1].
            %
            % Outputs:
            %   vol               - (double) Calculated volume of the polyhedral zone.
            %
            % Uses:
            %   tet_vol           - Calculates volume of a tetrahedron.
            %   cal_zone_ctr      - Calculates centroid of a 3D zone.
            %
            % Modifies:
            %   None.
            %
            % Original C declaration:
            % void cal_vol(int nface, int nnode, double *coords,
            %              int *nnode_for_face, int *nodelist_for_face,
            %              double *vol)
        
            % If polyhedron is a tetrahedron
            if nnode == 4
                vol = obj.tet_vol(coords(1, :), coords(2, :), coords(3, :), coords(4, :));
            else
                % Calculate the centroid of the zone
                ctr = obj.cal_zone_ctr(nface, nnode, coords, nnode_for_face, nodelist_for_face);
        
                % Initialize volume
                vol = 0.0;
                offset = 1;
        
                % Calculate volume by iterating over each face
                for f = 1:nface
                    nn = nnode_for_face(f);
                    nodes = nodelist_for_face(offset:offset + nn - 1);
        
                    % Calculate face centroid
                    fctr = mean(coords(nodes, :), 1);
        
                    % Calculate volume for each side formed by face centroid, zone centroid,
                    % and edge nodes on the face
                    for k = 1:nn
                        k1 = mod(k, nn) + 1;
                        n0 = nodes(k);
                        n1 = nodes(k1);
                        c0 = coords(n0, :);
                        c1 = coords(n1, :);
        
                        % Calculate the tetrahedral volume for the side formed by ctr, fctr, c0, and c1
                        dvol = obj.tet_vol(ctr, fctr, c0, c1);
                        vol = vol + dvol;
                    end
        
                    % Update offset for the next face
                    offset = offset + nn;
                end
            end
        end

        function [ifintersected_face, nodelist_intersected_face, nodes2_in_between] = ...
            interface_face_intersect(~, nface, nedge, nnode, nnode_for_face, ...
                                     nodelist_for_face, edgelist_for_face, ...
                                     nodelist_for_edge, facelist_for_edge, ...
                                     nedge_for_node, edgelist_for_node, ...
                                     nnode_tot, nnode_interface, ...
                                     nodelist_for_interface, interface_node_to_edge)
        
            % interface_face_intersect - Determine intersection of a 3D interface with each face.
            %
            % This function calculates the intersection points between a 3D interface 
            % and each face of a polyhedral zone. It identifies intersected faces and 
            % determines surrounding nodes for each interface node.
            %
            % Inputs:
            %   nface                    - (int) Number of faces in the polyhedral zone.
            %   nedge                    - (int) Number of edges in the polyhedral zone.
            %   nnode                    - (int) Number of nodes in the polyhedral zone.
            %   nnode_for_face           - (int array) Number of nodes per face, size: [nface, 1].
            %   nodelist_for_face        - (int array) Node indices for each face, 
            %                              size: [sum(nnode_for_face), 1].
            %   edgelist_for_face        - (int array) Edge indices for each face, 
            %                              size: [sum(nnode_for_face), 1].
            %   nodelist_for_edge        - (int array) Node pairs for each edge, size: [nedge, 2].
            %   facelist_for_edge        - (int array) Face indices associated with each edge, 
            %                              size: [2*nedge, 1].
            %   nedge_for_node           - (int array) Number of edges associated with each node, 
            %                              size: [nnode, 1].
            %   edgelist_for_node        - (int array) List of edges associated with each node, 
            %                              size: [sum(nedge_for_node), 1].
            %   nnode_tot                - (int) Total number of nodes, including interface nodes.
            %   nnode_interface          - (int) Number of nodes on the interface.
            %   nodelist_for_interface   - (int array) Node indices on the interface, 
            %                              size: [nnode_interface, 1].
            %   interface_node_to_edge   - (int array) Edge indices for each interface node, 
            %                              size: [nnode_interface, 1].
            %
            % Outputs:
            %   ifintersected_face       - (int array) Indicator for face intersection with 
            %                              the interface (1 if intersected, 0 otherwise), 
            %                              size: [nface, 1].
            %   nodelist_intersected_face - (int array) List of intersecting nodes for each face 
            %                               that intersects the interface, size: [2*nface, 1]. 
            %                               -1 if no intersection.
            %   nodes2_in_between        - (int array) Surrounding nodes for each interface node, 
            %                              size: [2*nnode_interface, 1].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   ifintersected_face         - Marks intersected faces with the interface.
            %   nodelist_intersected_face  - Stores the intersecting node indices for each face.
            %   nodes2_in_between          - Stores surrounding node indices for each interface node.
            %
            % Original C declaration:
            % void interface_face_intersect(int nface, int nedge, int nnode,
            %                               int *nnode_for_face, int *nodelist_for_face, int *edgelist_for_face,
            %                               int *nodelist_for_edge, int *facelist_for_edge,
            %                               int *nedge_for_node, int *edgelist_for_node, 
            %                               int nnode_tot, int nnode_interface, int *nodelist_for_interface, 
            %                               int *interface_node_to_edge, int *ifintersected_face,
            %                               int *nodelist_intersected_face, int **nodes2_in_between);
        
            % Initialize output variables
            ifintersected_face = zeros(nface, 1);
            nodelist_intersected_face = -1 * ones(2 * nface, 1);
            nodes2_in_between = zeros(2 * nnode_interface, 1);
        
            % Early exit if there are fewer than 3 interface nodes or if interface is a face
            if nnode_interface < 3 || nnode_tot == nnode
                return;
            end
        
            % Allocate edges_for_node to store edge lists for each node
            edges_for_node = cell(nnode, 1);
            offset = 1;
            for i = 1:nnode
                edges_for_node{i} = edgelist_for_node(offset:offset + nedge_for_node(i) - 1);
                offset = offset + nedge_for_node(i);
            end
        
            % Mark all faces as initially non-intersected
            nnode_intersecred_face = zeros(nface, 1);
        
            % Process each interface node
            for i = 1:nnode_interface
                n = nodelist_for_interface(i);
                e = interface_node_to_edge(i);
        
                if e < 0  % n is an original node
                    nodes2_in_between(2 * i - 1:2 * i) = n;
                    edges = edges_for_node{n};
                    
                    % Mark intersected faces
                    for k = 1:length(edges)
                        mye = edges(k);
                        faces = facelist_for_edge(2 * mye - 1:2 * mye);
        
                        for j = 1:2
                            f = faces(j);
                            nn = nnode_intersecred_face(f);
                            nodes = nodelist_intersected_face(2 * f - 1:2 * f);
        
                            if ~any(nodes(1:nn) == n)
                                assert(nn < 2);
                                nodelist_intersected_face(2 * f - 2 + nn + 1) = n;
                                nnode_intersecred_face(f) = nnode_intersecred_face(f) + 1;
                            end
                        end
                    end
                else  % n is an intersection
                    nodes2_in_between(2 * i - 1:2 * i) = nodelist_for_edge(2 * e - 1:2 * e);
                    faces = facelist_for_edge(2 * e - 1:2 * e);
        
                    % Mark intersected faces
                    for j = 1:2
                        f = faces(j);
                        nn = nnode_intersecred_face(f);
                        nodes = nodelist_intersected_face(2 * f - 1:2 * f);
        
                        if ~any(nodes(1:nn) == n)
                            assert(nn < 2);
                            nodelist_intersected_face(2 * f - 2 + nn + 1) = n;
                            nnode_intersecred_face(f) = nnode_intersecred_face(f) + 1;
                        end
                    end
                end
            end
        
            % Finalize face intersection information
            for i = 1:nface
                if nnode_intersecred_face(i) < 2
                    nnode_intersecred_face(i) = 0;
                    nodelist_intersected_face(2 * i - 1) = -1;
                else
                    nnode_intersecred_face(i) = 1;
                end
            end
        end

        function [xmin, dxmax] = get_scaling(~, nnode, coords)
            % get_scaling - Calculate scaling factors for a 3D mesh based on node coordinates.
            %
            % This function calculates the minimum and maximum extents of the coordinates 
            % of a 3D mesh along each axis, and determines the largest span among all axes.
            %
            % Inputs:
            %   nnode   - (int) Number of nodes in the mesh.
            %   coords  - (double array) Coordinates of all nodes, size: [nnode, 3].
            %
            % Outputs:
            %   xmin    - (double array) Minimum coordinates along each axis, size: [1, 3].
            %   dxmax   - (double) Maximum span across all axes.
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   None.
            %
            % Original C declaration:
            % void get_scaling(int nnode, double *coords, double *xmin, double *dxmax)
        
            % Ensure at least one node is present
            assert(nnode > 0);
        
            % Initialize minimum and maximum coordinates
            xmin = coords(1, :);
            xmax = coords(1, :);
        
            % Iterate over all nodes to find min and max values for each axis
            for n = 2:nnode
                c = coords(n, :);
                xmin = min(xmin, c);
                xmax = max(xmax, c);
            end
        
            % Calculate the maximum span across all dimensions
            dxmax = max(xmax - xmin);
        end

        function norm_ea_face = cal_norm_of_face(obj, nface, nnode, coords, nnode_ea_face, ...
                                                 nodelist_for_face, edgelist_for_face, ...
                                                 interior_point)
            % cal_norm_of_face - Calculate outward normals for each face of a 3D polyhedron.
            %
            % This function calculates the outward-facing normal vector for each face of 
            % a polyhedron defined by nodes and faces. If `interior_point` is not provided, 
            % the function calculates the centroid of the polyhedron.
            %
            % Inputs:
            %   obj                - (c_vof3d) Encapsulating class object
            %   nface              - (int) Number of faces in the polyhedral zone.
            %   nnode              - (int) Number of nodes in the polyhedral zone.
            %   coords             - (double array) Coordinates of all nodes, size: [nnode, 3].
            %   nnode_ea_face      - (int array) Number of nodes per face, size: [nface, 1].
            %   nodelist_for_face  - (int array) Node indices for each face, 
            %                        size: [sum(nnode_ea_face), 1].
            %   edgelist_for_face  - (int array) Edge indices for each face, 
            %                        size: [sum(nnode_ea_face), 1].
            %   interior_point     - (double array) Coordinates of an interior point, 
            %                        size: [1, 3]. If empty, the centroid is computed.
            %
            % Outputs:
            %   norm_ea_face       - (double array) Outward normals for each face, 
            %                        size: [nface, 3].
            %
            % Uses:
            %   get_interior_point - Computes the centroid if `interior_point` is not provided.
            %
            % Modifies:
            %   norm_ea_face       - Stores the outward normals for each face.
            %
            % Original C declaration:
            % void cal_norm_of_face(int nface, int nnode, double *coords,
            %                       int *nnode_ea_face, int *nodelist_for_face,
            %                       int *edgelist_for_face, double *interior_point,
            %                       double *norm_ea_face)
        
            % Initialize variables
            dim = 3;
            small = 1e-10;
            norm_ea_face = zeros(nface, dim);
        
            % Calculate the centroid if interior_point is not provided
            if isempty(interior_point)
                interior_point = obj.get_interior_point(dim, nface, nnode, coords, ...
                                                    nnode_ea_face, nodelist_for_face);
            end
        
            % Process each face to calculate the normal
            offset = 1;
            for f = 1:nface
                nn = nnode_ea_face(f);
                nodes = nodelist_for_face(offset:offset + nn - 1);
                
                % Calculate face centroid
                ctrf = mean(coords(nodes, :), 1);
        
                % Calculate the face normal
                found = false;
                for i = 1:nn
                    i1 = mod(i, nn) + 1;
                    n = nodes(i);
                    n1 = nodes(i1);
                    r0 = coords(n, :) - ctrf;
                    r1 = coords(n1, :) - ctrf;
        
                    % Cross product for the normal
                    norm = cross(r0, r1);
                    if norm * norm' > small
                        norm = norm / norm(norm);  % Normalize
                        found = true;
                        break;
                    end
                end
        
                % Default normal if cross product is zero
                if ~found
                    norm = [1.0, 0.0, 0.0];
                end
        
                % Ensure outward orientation based on interior_point
                n = nodes(1);
                if dot(coords(n, :) - interior_point, norm) < 0
                    norm = -norm;
                    
                    % Reverse order of nodes and edges for outward-facing normal
                    edges = edgelist_for_face(offset:offset + nn - 1);
                    nodes = flip(nodes);
                    edges = flip(edges);
                    nodelist_for_face(offset:offset + nn - 1) = nodes;
                    edgelist_for_face(offset:offset + nn - 1) = edges;
                end
        
                % Assign normal to the output array
                norm_ea_face(f, :) = norm;
                
                % Move offset for the next face
                offset = offset + nn;
            end
        end

        function [nnode_tot, coords_tot, nnode_interface, nodelist_for_interface, interfacenode_to_edge] = ...
            find_interface_z(obj, nface, nedge, nnode, coords, nnode_ea_face, nodelist_for_face, ...
                             edgelist_for_face, nodelist_for_edge, nedge_ea_node, ...
                             edgelist_for_node, facelist_for_edge, node_loc, ds)
            % find_interface_z - Locate a material interface in a 3D polyhedral zone at a specific z-level.
            %
            % This function identifies a material interface in a 3D polyhedral zone, 
            % based on a given z-coordinate (`ds`). It calculates the intersection of 
            % the interface with edges and nodes in the polyhedron.
            %
            % Inputs:
            %   obj                   - (c_vof3d) Encapsulation class object
            %   nface                 - (int) Number of faces in the polyhedral zone.
            %   nedge                 - (int) Number of edges in the polyhedral zone.
            %   nnode                 - (int) Number of nodes in the polyhedral zone.
            %   coords                - (double array) Coordinates of all nodes, size: [nnode, 3].
            %   nnode_ea_face         - (int array) Number of nodes per face, size: [nface, 1].
            %   nodelist_for_face     - (int array) Node indices for each face, size: [sum(nnode_ea_face), 1].
            %   edgelist_for_face     - (int array) Edge indices for each face, size: [sum(nnode_ea_face), 1].
            %   nodelist_for_edge     - (int array) Node pairs for each edge, size: [nedge, 2].
            %   nedge_ea_node         - (int array) Number of edges per node, size: [nnode, 1].
            %   edgelist_for_node     - (int array) List of edges associated with each node, 
            %                           size: [sum(nedge_ea_node), 1].
            %   facelist_for_edge     - (int array) Face indices associated with each edge, 
            %                           size: [2*nedge, 1].
            %   node_loc              - (int array) Node classification based on z-position 
            %                           relative to the interface (`ds`). Values: 
            %                           -1 (below), 0 (on interface), 1 (above), size: [nnode, 1].
            %   ds                    - (double) Z-coordinate of the material interface.
            %
            % Outputs:
            %   nnode_tot             - (int) Total number of nodes, including interface nodes.
            %   coords_tot            - (double array) Coordinates of all nodes after interface 
            %                           inclusion, size: [nnode_tot, 3].
            %   nnode_interface       - (int) Number of nodes on the material interface.
            %   nodelist_for_interface - (int array) Node indices on the material interface, 
            %                            size: [nnode_interface, 1].
            %   interfacenode_to_edge - (int array) Edge indices associated with each interface 
            %                           node, size: [nnode_interface, 1].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   nnode_tot              - Updated to reflect the total node count after adding 
            %                            interface nodes.
            %   coords_tot             - Extended to include coordinates of interface nodes.
            %   nnode_interface        - Updated with the number of nodes on the interface.
            %   nodelist_for_interface - Updated to include node indices on the interface.
            %   interfacenode_to_edge  - Stores edges associated with each interface node.
            %
            % Original C declaration:
            % void find_interface_z(int nface, int nedge, int nnode, double *coords,
            %                       int *nnode_ea_face, int *nodelist_for_face,
            %                       int *edgelist_for_face, int *nodelist_for_edge,
            %                       int *nedge_ea_node, int *edgelist_for_node,
            %                       int *facelist_for_edge, int *node_loc, double ds,
            %                       int *nnode_tot, double **coords_tot,
            %                       int *nnode_interface, int *nodelist_for_interface,
            %                       int *interfacenode_to_edge);
        
            % Constants and initialization
            dim = 3;
            small = 1e-10;  % Small tolerance for distinguishing points
            dzmin = 1e-8;   % Minimum z-distance for edge detection
            
            % Initialize outputs and temporary storage
            nnode_interface = 0;
            nnode_tot = nnode;
            edge_taken = zeros(nedge, 1);
            coords_new = zeros(nedge, dim);
            nodelist_for_interface = zeros(nedge, 1);
            interfacenode_to_edge = zeros(nedge, 1);
        
            % Initialize edge list for each node
            edges_ea_node = cell(nnode, 1);
            offset = 1;
            for i = 1:nnode
                edges_ea_node{i} = edgelist_for_node(offset:offset + nedge_ea_node(i) - 1);
                offset = offset + nedge_ea_node(i);
            end
        
            % Mark edges as taken if they are close in z-coordinate
            nodes = nodelist_for_edge;
            for e = 1:nedge
                n0 = nodes(2 * e - 1);
                n1 = nodes(2 * e);
                c0 = coords(n0, :);
                c1 = coords(n1, :);
                
                if abs(c0(3) - c1(3)) < dzmin
                    edge_taken(e) = 1;
                end
            end
        
            % Process interface nodes
            for i = 1:nnode_interface
                n = nodelist_for_interface(i);
                node_loc(n) = 0;
        
                % Mark all edges connected to this node as taken
                edges = edges_ea_node{n};
                edge_taken(edges) = 1;
                interfacenode_to_edge(i) = -1;
            end
        
            % Identify new nodes created at interface edges
            ncoord_new = 0;
            for e = 1:nedge
                if edge_taken(e), continue; end
        
                n0 = nodelist_for_edge(2 * e - 1);
                n1 = nodelist_for_edge(2 * e);
                if node_loc(n0) * node_loc(n1) == 0
                    edge_taken(edges_ea_node{n0}) = 1;
                    edge_taken(edges_ea_node{n1}) = 1;
                end
        
                if node_loc(n0) * node_loc(n1) > 0, continue; end
        
                % Calculate intersection point at ds level
                c0 = coords(n0, :);
                c1 = coords(n1, :);
                point = obj.line_intersect_z(c0, c1, ds);
                coords_new(ncoord_new + 1, :) = point;
        
                % Update interface node information
                nodelist_for_interface(nnode_interface + 1) = nnode + ncoord_new;
                interfacenode_to_edge(nnode_interface + 1) = e;
                nnode_interface = nnode_interface + 1;
                ncoord_new = ncoord_new + 1;
        
                edge_taken(e) = 1;
            end
        
            % Update total nodes and concatenate coordinates with interface nodes
            nnode_tot = nnode + ncoord_new;
            coords_tot = [coords; coords_new(1:ncoord_new, :)];
        
            % Sort the interface nodes if more than three nodes are present
            if nnode_interface > 3
                coords_work = zeros(nnode_interface, 2);
                node_new2old = (1:nnode_interface)';
        
                for i = 1:nnode_interface
                    n0 = nodelist_for_interface(i);
                    if n0 <= nnode
                        coords_work(i, :) = coords(n0, 1:2);
                    else
                        coords_work(i, :) = coords_new(n0 - nnode, 1:2);
                    end
                end
        
                [nodelist_for_interface, interfacenode_to_edge] = obj.sort_nodes( ...
                    nnode_interface, coords_work, interfacenode_to_edge);
        
                % Map sorted nodes back to original indices
                nodelist_for_interface = node_new2old(nodelist_for_interface);
            end
        end

        function [nnode_tot, coords_tot, nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, ...
          face2old_lower, nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, ...
          face2old_upper, faceindexlow_of_upper_interface, faceindexhgh_of_upper_interface] = ...
          divide_z(obj, nface, nnode, coords, ds, node_loc, nnode_ea_face, nodelist_for_face, ...
                   nedge, edgelist_for_face, nodelist_for_edge, nedge_ea_node, edgelist_for_node, ...
                   facelist_for_edge, norm_ea_face, nnode_input, coords_input, nnode_interface, ...
                   nodes_interface, interface_node_to_edge)
            % divide_z - Divide a polyhedron into lower and upper parts along a z-plane.
            %
            % This function divides a 3D polyhedron into lower and upper polyhedrons based 
            % on a given z-coordinate (`ds`). It identifies and organizes the faces, edges, 
            % and nodes for each resulting polyhedron, including the interface.
            %
            % Inputs:
            %   nface                             - (int) Number of faces in the polyhedral zone.
            %   nnode                             - (int) Number of nodes in the polyhedral zone.
            %   coords                            - (double array) Coordinates of all nodes, 
            %                                       size: [nnode, 3].
            %   ds                                - (double) Z-coordinate of the division plane.
            %   node_loc                          - (int array) Node classification based on 
            %                                       z-position relative to `ds`.
            %                                       Values: -1 (below), 0 (on interface), 
            %                                       1 (above), size: [nnode, 1].
            %   nnode_ea_face                     - (int array) Number of nodes per face, 
            %                                       size: [nface, 1].
            %   nodelist_for_face                 - (int array) Node indices for each face, 
            %                                       size: [sum(nnode_ea_face), 1].
            %   nedge                             - (int) Number of edges in the polyhedral zone.
            %   edgelist_for_face                 - (int array) Edge indices for each face, 
            %                                       size: [sum(nnode_ea_face), 1].
            %   nodelist_for_edge                 - (int array) Node pairs for each edge, 
            %                                       size: [nedge, 2].
            %   nedge_ea_node                     - (int array) Number of edges per node, 
            %                                       size: [nnode, 1].
            %   edgelist_for_node                 - (int array) List of edges associated with 
            %                                       each node, size: [sum(nedge_ea_node), 1].
            %   facelist_for_edge                 - (int array) Face indices associated with 
            %                                       each edge, size: [2*nedge, 1].
            %   norm_ea_face                      - (double array) Outward normals for each face, 
            %                                       size: [nface, 3].
            %   nnode_input                       - (int) Number of input nodes.
            %   coords_input                      - (double array) Input node coordinates, 
            %                                       size: [nnode_input, 3].
            %   nnode_interface                   - (int) Number of nodes on the interface.
            %   nodes_interface                   - (int array) Node indices on the interface, 
            %                                       size: [nnode_interface, 1].
            %   interface_node_to_edge            - (int array) Edge indices for each interface 
            %                                       node, size: [nnode_interface, 1].
            %   faceindexlow_of_upper_interface   - (int) Index of the interface face in the lower 
            %                                       polyhedron.
            %   faceindexhgh_of_upper_interface   - (int) Index of the interface face in the upper 
            %                                       polyhedron.
            %
            % Outputs:
            %   nnode_tot                         - (int) Total number of nodes, including 
            %                                       interface nodes.
            %   coords_tot                        - (double array) Coordinates of all nodes 
            %                                       after division, size: [nnode_tot, 3].
            %   nface_lower                       - (int) Number of faces in the lower polyhedron.
            %   nnode_ea_face_lower               - (int array) Number of nodes per face for 
            %                                       the lower polyhedron, size: [nface_lower, 1].
            %   nodelist_for_face_lower           - (int array) Node indices for each face in 
            %                                       the lower polyhedron, 
            %                                       size: [sum(nnode_ea_face_lower), 1].
            %   face2old_lower                    - (int array) Mapping from lower polyhedron 
            %                                       faces to original mesh faces, size: [nface_lower, 1].
            %   nface_upper                       - (int) Number of faces in the upper polyhedron.
            %   nnode_ea_face_upper               - (int array) Number of nodes per face for 
            %                                       the upper polyhedron, size: [nface_upper, 1].
            %   nodelist_for_face_upper           - (int array) Node indices for each face in 
            %                                       the upper polyhedron, 
            %                                       size: [sum(nnode_ea_face_upper), 1].
            %   face2old_upper                    - (int array) Mapping from upper polyhedron 
            %                                       faces to original mesh faces, size: [nface_upper, 1].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   nnode_tot, coords_tot             - Updated with the total node count and 
            %                                       coordinates after division.
            %   nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower - 
            %                                       Define the lower polyhedron.
            %   nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper - 
            %                                       Define the upper polyhedron.
            %
            % Original C declaration:
            % void divide_z(int nface, int nnode, double *coords,
            %               double ds, int *node_loc,
            %               int *nnode_ea_face, int *nodelist_for_face,
            %               int nedge, int *edgelist_for_face, int *nodelist_for_edge,
            %               int *nedge_ea_node, int *edgelist_for_node,
            %               int *facelist_for_edge, double **norm_ea_face,
            %               int nnode_input, double *coords_input,
            %               int nnode_interface, int *nodes_interface,
            %               int *interface_node_to_edge,
            %               int *faceindexlow_of_upper_interface, int *faceindexhgh_of_upper_interface,
            %               int *nnode_tot, double **coords_tot,
            %               int *nface_lower, int *nnode_ea_face_lower, 
            %               int *nodelist_for_face_lower, int *face2old_lower,
            %               int *nface_upper, int *nnode_ea_face_upper, 
            %               int *nodelist_for_face_upper, int *face2old_upper)


            % Constants and initialization
            dim = 3;
            nnode_new = 0;
            mxnnode_new = 2 * nnode + 1;
            coords_new = zeros(mxnnode_new, dim);
            edge_loc = zeros(nedge, 1);
            edge_to_node = -1 * ones(nedge, 1);
            node_taken = zeros(nnode_interface, 1);
            faceindexlow_of_upper_interface = -1;
            faceindexhgh_of_upper_interface = -1;
        
            % Process interface nodes
            for i = 1:nnode_interface
                n0 = nodes_interface(i);
                if n0 <= nnode
                    continue;
                end
        
                e = interface_node_to_edge(i);
                edge_loc(e) = 2;  % Mark edge as crossing the interface
                nodes = nodelist_for_edge(2 * e - 1:2 * e);
                n0 = nodes(1);
                n1 = nodes(2);
                c0 = coords(n0, :);
                c1 = coords(n1, :);
        
                % Calculate intersection point
                point = line_intersect_z(c0, c1, ds);
                coords_new(nnode_new + 1, :) = point;
                edge_to_node(e) = nnode + nnode_new;
                nodes_interface(i) = nnode + nnode_new;
                nnode_new = nnode_new + 1;
        
                % Adjust node_loc
                if node_loc(n0) == 0
                    node_loc(n0) = -1 + 2 * (node_loc(n1) > 0);
                elseif node_loc(n1) == 0
                    node_loc(n1) = -1 + 2 * (node_loc(n0) > 0);
                end
            end
        
            % Process edges to classify them
            for e = 1:nedge
                if edge_loc(e) == 2
                    continue;
                end
                nodes = nodelist_for_edge(2 * e - 1:2 * e);
                n0 = nodes(1);
                n1 = nodes(2);
        
                if node_loc(n0) == 0 && node_loc(n1) == 0
                    edge_loc(e) = 0;  % On the interface
                elseif node_loc(n0) >= 0 && node_loc(n1) >= 0
                    edge_loc(e) = 1;  % Above the interface
                elseif node_loc(n0) <= 0 && node_loc(n1) <= 0
                    edge_loc(e) = -1; % Below the interface
                end
            end
        
            % Check if the interface is an original face
            interface_is_face = 0;
            if nnode_new == 0
                offset = 1;
                for f = 1:nface
                    nn = nnode_ea_face(f);
                    if nn ~= nnode_interface
                        offset = offset + nn;
                        continue;
                    end
        
                    nodes = nodelist_for_face(offset:offset + nn - 1);
                    node_taken(1:nn) = 0;
                    found = true;
        
                    for i = 1:nnode_interface
                        n = nodes_interface(i);
                        matched = false;
        
                        for k = 1:nn
                            if ~node_taken(k) && nodes(k) == n
                                node_taken(k) = 1;
                                matched = true;
                                break;
                            end
                        end
        
                        if ~matched
                            found = false;
                            break;
                        end
                    end
        
                    if found
                        interface_is_face = f;  % 1-based indexing for MATLAB
                        break;
                    end
        
                    offset = offset + nn;
                end
            end
        
            % Generate lower and upper polyhedrons
            [nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower, ...
             nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper, ...
             faceindexlow_of_upper_interface, faceindexhgh_of_upper_interface] = ...
                obj.get_lower_upper_polyhedron(nface, nedge, nnode, nnode_ea_face, ...
                                           nodelist_for_face, edgelist_for_face, ...
                                           nodelist_for_edge, nnode_new, edge_to_node, ...
                                           nnode_interface, nodes_interface, node_loc, ...
                                           edge_loc, interface_is_face);
        
            % Update total nodes and concatenate coordinates
            nnode_tot = nnode + nnode_new;
            coords_tot = [coords; coords_new(1:nnode_new, :)];
        
            % Validate polyhedron faces
            if nface_upper < 4
                nface_upper = 0;
            end
            if nface_lower < 4
                nface_lower = 0;
            end
        end

        function [ds, nnode_tot, coords_tot, nnode_for_interface, nodelist_for_interface, ...
                  interface_node_to_edge, faceindexlow_of_upper_interface, ...
                  faceindexhgh_of_upper_interface, nface_lower, nnode_ea_face_lower, ...
                  nodelist_for_face_lower, face2old_lower, nface_upper, ...
                  nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper] = ...
                  bisection_z(obj, vf_to_match, volume, ds_lower, ds_upper, vf_lower, ...
                              vf_upper, nface, nedge, nnode, coords, nnode_ea_face, ...
                              nodelist_for_face, edgelist_for_face, nodelist_for_edge)
            % bisection_z - Find the z-coordinate that matches a target volume fraction.
            %
            % This function determines the z-coordinate (`ds`) at which a material interface 
            % divides a polyhedron such that the volume occupied by one material matches a 
            % given target volume fraction (`vf_to_match`).
            %
            % Inputs:
            %   vf_to_match                     - (double) Target volume fraction to match.
            %   volume                          - (double) Total volume of the polyhedron.
            %   ds_lower                        - (double) Lower bound for the z-coordinate.
            %   ds_upper                        - (double) Upper bound for the z-coordinate.
            %   vf_lower                        - (double) Volume fraction at `ds_lower`.
            %   vf_upper                        - (double) Volume fraction at `ds_upper`.
            %   nface                           - (int) Number of faces in the polyhedral zone.
            %   nedge                           - (int) Number of edges in the polyhedral zone.
            %   nnode                           - (int) Number of nodes in the polyhedral zone.
            %   coords                          - (double array) Coordinates of all nodes, 
            %                                     size: [nnode, 3].
            %   nnode_ea_face                   - (int array) Number of nodes per face, 
            %                                     size: [nface, 1].
            %   nodelist_for_face               - (int array) Node indices for each face, 
            %                                     size: [sum(nnode_ea_face), 1].
            %   edgelist_for_face               - (int array) Edge indices for each face, 
            %                                     size: [sum(nnode_ea_face), 1].
            %   nodelist_for_edge               - (int array) Node pairs for each edge, 
            %                                     size: [nedge, 2].
            %
            % Outputs:
            %   ds                              - (double) Computed z-coordinate that matches 
            %                                     the target volume fraction.
            %   nnode_tot                       - (int) Total number of nodes, including 
            %                                     interface nodes.
            %   coords_tot                      - (double array) Coordinates of all nodes after 
            %                                     adding interface nodes, size: [nnode_tot, 3].
            %   nnode_for_interface             - (int) Number of nodes on the material interface.
            %   nodelist_for_interface          - (int array) Node indices on the material 
            %                                     interface, size: [nnode_for_interface, 1].
            %   interface_node_to_edge          - (int array) Edge indices for each interface 
            %                                     node, size: [nnode_for_interface, 1].
            %   faceindexlow_of_upper_interface - (int) Index of the interface face in the lower 
            %                                     polyhedron.
            %   faceindexhgh_of_upper_interface - (int) Index of the interface face in the upper 
            %                                     polyhedron.
            %   nface_lower                     - (int) Number of faces in the lower polyhedron.
            %   nnode_ea_face_lower             - (int array) Number of nodes per face in the 
            %                                     lower polyhedron, size: [nface_lower, 1].
            %   nodelist_for_face_lower         - (int array) Node indices for each face in the 
            %                                     lower polyhedron, 
            %                                     size: [sum(nnode_ea_face_lower), 1].
            %   face2old_lower                  - (int array) Mapping from lower polyhedron 
            %                                     faces to original mesh faces, size: [nface_lower, 1].
            %   nface_upper                     - (int) Number of faces in the upper polyhedron.
            %   nnode_ea_face_upper             - (int array) Number of nodes per face in the 
            %                                     upper polyhedron, size: [nface_upper, 1].
            %   nodelist_for_face_upper         - (int array) Node indices for each face in the 
            %                                     upper polyhedron, 
            %                                     size: [sum(nnode_ea_face_upper), 1].
            %   face2old_upper                  - (int array) Mapping from upper polyhedron 
            %                                     faces to original mesh faces, size: [nface_upper, 1].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   ds                              - Updated to the computed z-coordinate.
            %   nnode_tot, coords_tot           - Updated with total node count and coordinates 
            %                                     after division.
            %   nnode_for_interface, nodelist_for_interface, interface_node_to_edge - Updated 
            %                                     with interface node and edge information.
            %   nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower - 
            %                                     Define the lower polyhedron.
            %   nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper - 
            %                                     Define the upper polyhedron.
            %
            % Original C declaration:
            % void bisection_z(double vf_to_match, double volume,
            %                  double ds_lower, double ds_upper,
            %                  double vf_lower, double vf_upper,
            %                  int nface, int nedge, int nnode,
            %                  double *coords,
            %                  int *nnode_ea_face, int *nodelist_for_face,
            %                  int *edgelist_for_face, int *nodelist_for_edge,
            %                  double *ds, int *nnode_tot, double **coords_tot,
            %                  int *nnode_for_interface, int **nodelist_for_interface, 
            %                  int **interface_node_to_edge,
            %                  int *faceindexlow_of_upper_interface, 
            %                  int *faceindexhgh_of_upper_interface,
            %                  int *nface_lower, int *nnode_ea_face_lower, 
            %                  int *nodelist_for_face_lower, int *face2old_lower,
            %                  int *nface_upper, int *nnode_ea_face_upper, 
            %                  int *nodelist_for_face_upper, int *face2old_upper);



            % Constants and initialization
            dim = 3;
            niter_mx = 50;  % Maximum iterations
            accuracy = 1e-6; % Convergence threshold
            dzmin = 1e-8;   % Minimum tolerance for z-level comparison
            ds0 = ds_lower;
            ds1 = ds_upper;
            vf0 = vf_lower;
            vf1 = vf_upper;
            vf = vf0;
            ds = 0.5 * (ds0 + ds1);
            err = 1.0;
            niter = 0;
            done = false;
        
            % Allocate memory for temporary storage
            node_loc = zeros(nnode, 1);
            coords_tot = [];
            facelist_for_edge = zeros(2 * nedge, 1);
            [nedge_ea_node, edgelist_for_node] = get_faces_ea_edge(nface, nedge, ...
                nnode, nnode_ea_face, nodelist_for_face, edgelist_for_face, ...
                nodelist_for_edge, facelist_for_edge);
        
            norm_ea_face = obj.cal_norm_of_face(nface, nnode, coords, ...
                                            nnode_ea_face, nodelist_for_face, ...
                                            edgelist_for_face, []);
        
            % Interface variables
            nnode_for_interface_previous = 0;
            nodelist_for_interface_previous = [];
            nodelist_for_interface = zeros(nedge + nnode, 1);
            interface_node_to_edge = zeros(nedge + nnode, 1);
        
            % Begin bisection loop
            while ~done && (niter < niter_mx)
                % Classify nodes relative to the current z-level
                for i = 1:nnode
                    if abs(coords(i, 3) - ds) <= dzmin
                        node_loc(i) = 0; % On interface
                    elseif coords(i, 3) > ds
                        node_loc(i) = 1; % Above interface
                    else
                        node_loc(i) = -1; % Below interface
                    end
                end
        
                % Find interface nodes
                [nnode_all, coords_all, nnode_for_interface, ...
                    nodelist_for_interface, interface_node_to_edge] = ...
                    obj.find_interface_z(nface, nedge, nnode, coords, nnode_ea_face, ...
                                     nodelist_for_face, edgelist_for_face, ...
                                     nodelist_for_edge, nedge_ea_node, ...
                                     edgelist_for_node, facelist_for_edge, ...
                                     node_loc, ds);
        
                if nnode_for_interface < 3
                    % Handle degenerate interface cases
                    if niter == 0
                        nnode_tot = nnode;
                        coords_tot = coords;
                        sz = sum(nnode_ea_face);
        
                        if vf_to_match < 0.5
                            nface_lower = 0;
                            nface_upper = nface;
                            face2old_upper = 1:nface;
                            nnode_ea_face_upper = nnode_ea_face;
                            nodelist_for_face_upper = nodelist_for_face;
                        else
                            nface_upper = 0;
                            nface_lower = nface;
                            face2old_lower = 1:nface;
                            nnode_ea_face_lower = nnode_ea_face;
                            nodelist_for_face_lower = nodelist_for_face;
                        end
                    else
                        % Revert to previous values
                        ds = ds_previous;
                        nnode_for_interface = nnode_for_interface_previous;
                        nodelist_for_interface = nodelist_for_interface_previous;
                    end
                    done = true;
                else
                    % Divide polyhedron and compute volume
                    [nnode_tot, coords_tot, nface_lower, nnode_ea_face_lower, ...
                        nodelist_for_face_lower, face2old_lower, nface_upper, ...
                        nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper] = ...
                        obj.divide_z(nface, nnode, coords, ds, node_loc, nnode_ea_face, ...
                                 nodelist_for_face, nedge, edgelist_for_face, ...
                                 nodelist_for_edge, nedge_ea_node, edgelist_for_node, ...
                                 facelist_for_edge, norm_ea_face, nnode_all, coords_all, ...
                                 nnode_for_interface, nodelist_for_interface, ...
                                 interface_node_to_edge);
        
                    % Calculate volume of the lower polyhedron
                    vol = obj.cal_vol(nface_lower, nnode_tot, coords_tot, ...
                                  nnode_ea_face_lower, nodelist_for_face_lower);
                    vf = vol / volume;
        
                    % Check convergence
                    err = abs(vf - vf_to_match);
                    if err <= obj.accuracy || ...
                       (ds1 - ds0) / (ds_upper - ds_lower) <= obj.accuracy
                        done = true;
                    else
                        if vf > vf_to_match
                            vf1 = vf;
                            ds1 = ds;
                        else
                            vf0 = vf;
                            ds0 = ds;
                        end
                        ds_previous = ds;
                        ds = 0.5 * (ds0 + ds1);
                        niter = niter + 1;
                    end
                    nnode_for_interface_previous = nnode_for_interface;
                    nodelist_for_interface_previous = nodelist_for_interface;
                end
        
                if exist('coords_all', 'var') && ~isempty(coords_all)
                    clear coords_all;
                end
            end
        
            % Free temporary resources
            clear facelist_for_edge nedge_ea_node edgelist_for_node norm_ea_face;
        end

        function [ds_lower, ds_upper, vf_lower, vf_upper, nnode_tot, coords_tot, ...
                  vol_matched, nnode_upper_interface, nodelist_upper_interface, ...
                  faceindexlow_of_upper_interface, faceindexhgh_of_upper_interface, ...
                  ifintersected_face, nodelist_intersected_face, nodes2_in_between, ...
                  nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower, ...
                  nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper] = ...
                  bounds_z(obj, vf_to_match, volume, nface, nnode, coords, nnode_ea_face, ...
                           nodelist_for_face, node_order_for_ds, ds_ea_node, my_nedge, ...
                           my_edgelist_for_face, my_nodelist_for_edge, my_nedge_ea_node, ...
                           my_edgelist_for_node, my_facelist_for_edge, my_norm_ea_face)
                
            % bisection_z - Find the z-coordinate that matches a target volume fraction.
            %
            % This function determines the z-coordinate (`ds`) at which a material interface 
            % divides a polyhedron such that the volume occupied by one material matches a 
            % given target volume fraction (`vf_to_match`).
            %
            % Inputs:
            %   vf_to_match                     - (double) Target volume fraction to match.
            %   volume                          - (double) Total volume of the polyhedron.
            %   ds_lower                        - (double) Lower bound for the z-coordinate.
            %   ds_upper                        - (double) Upper bound for the z-coordinate.
            %   vf_lower                        - (double) Volume fraction at `ds_lower`.
            %   vf_upper                        - (double) Volume fraction at `ds_upper`.
            %   nface                           - (int) Number of faces in the polyhedral zone.
            %   nedge                           - (int) Number of edges in the polyhedral zone.
            %   nnode                           - (int) Number of nodes in the polyhedral zone.
            %   coords                          - (double array) Coordinates of all nodes, 
            %                                     size: [nnode, 3].
            %   nnode_ea_face                   - (int array) Number of nodes per face, 
            %                                     size: [nface, 1].
            %   nodelist_for_face               - (int array) Node indices for each face, 
            %                                     size: [sum(nnode_ea_face), 1].
            %   edgelist_for_face               - (int array) Edge indices for each face, 
            %                                     size: [sum(nnode_ea_face), 1].
            %   nodelist_for_edge               - (int array) Node pairs for each edge, 
            %                                     size: [nedge, 2].
            %
            % Outputs:
            %   ds                              - (double) Computed z-coordinate that matches 
            %                                     the target volume fraction.
            %   nnode_tot                       - (int) Total number of nodes, including 
            %                                     interface nodes.
            %   coords_tot                      - (double array) Coordinates of all nodes after 
            %                                     adding interface nodes, size: [nnode_tot, 3].
            %   nnode_for_interface             - (int) Number of nodes on the material interface.
            %   nodelist_for_interface          - (int array) Node indices on the material 
            %                                     interface, size: [nnode_for_interface, 1].
            %   interface_node_to_edge          - (int array) Edge indices for each interface 
            %                                     node, size: [nnode_for_interface, 1].
            %   faceindexlow_of_upper_interface - (int) Index of the interface face in the lower 
            %                                     polyhedron.
            %   faceindexhgh_of_upper_interface - (int) Index of the interface face in the upper 
            %                                     polyhedron.
            %   nface_lower                     - (int) Number of faces in the lower polyhedron.
            %   nnode_ea_face_lower             - (int array) Number of nodes per face in the 
            %                                     lower polyhedron, size: [nface_lower, 1].
            %   nodelist_for_face_lower         - (int array) Node indices for each face in the 
            %                                     lower polyhedron, 
            %                                     size: [sum(nnode_ea_face_lower), 1].
            %   face2old_lower                  - (int array) Mapping from lower polyhedron 
            %                                     faces to original mesh faces, size: [nface_lower, 1].
            %   nface_upper                     - (int) Number of faces in the upper polyhedron.
            %   nnode_ea_face_upper             - (int array) Number of nodes per face in the 
            %                                     upper polyhedron, size: [nface_upper, 1].
            %   nodelist_for_face_upper         - (int array) Node indices for each face in the 
            %                                     upper polyhedron, 
            %                                     size: [sum(nnode_ea_face_upper), 1].
            %   face2old_upper                  - (int array) Mapping from upper polyhedron 
            %                                     faces to original mesh faces, size: [nface_upper, 1].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   ds                              - Updated to the computed z-coordinate.
            %   nnode_tot, coords_tot           - Updated with total node count and coordinates 
            %                                     after division.
            %   nnode_for_interface, nodelist_for_interface, interface_node_to_edge - Updated 
            %                                     with interface node and edge information.
            %   nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower - 
            %                                     Define the lower polyhedron.
            %   nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper - 
            %                                     Define the upper polyhedron.
            %
            % Original C declaration:
            % void bisection_z(double vf_to_match, double volume,
            %                  double ds_lower, double ds_upper,
            %                  double vf_lower, double vf_upper,
            %                  int nface, int nedge, int nnode,
            %                  double *coords,
            %                  int *nnode_ea_face, int *nodelist_for_face,
            %                  int *edgelist_for_face, int *nodelist_for_edge,
            %                  double *ds, int *nnode_tot, double **coords_tot,
            %                  int *nnode_for_interface, int **nodelist_for_interface, 
            %                  int **interface_node_to_edge,
            %                  int *faceindexlow_of_upper_interface, 
            %                  int *faceindexhgh_of_upper_interface,
            %                  int *nface_lower, int *nnode_ea_face_lower, 
            %                  int *nodelist_for_face_lower, int *face2old_lower,
            %                  int *nface_upper, int *nnode_ea_face_upper, 
            %                  int *nodelist_for_face_upper, int *face2old_upper);
        
            % Constants and initialization
            dim = 3;
            dzmin = 1e-8; % Minimum z-difference for comparison
            accuracy = 1e-6; % Convergence threshold
            vol_matched = 0;
            coords_tot = [];
            ds_ordered = zeros(nnode, 1);
            node_loc = zeros(nnode, 1);
            vf = 0.0;
            vf_lower = 0.0;
            vf_upper = 0.0;
        
            % Edge and face data
            if my_nedge <= 0
                % Compute edges and normals dynamically
                [nedge, edgelist_for_face, nodelist_for_edge] = ...
                    obj.get_edges(nface, nnode, nnode_ea_face, nodelist_for_face);
        
                norm_ea_face = obj.cal_norm_of_face(nface, nnode, coords, nnode_ea_face, ...
                                                nodelist_for_face, edgelist_for_face);
        
                [facelist_for_edge, nedge_ea_node, edgelist_for_node] = ...
                    obj.get_faces_ea_edge(nface, nedge, nnode, nnode_ea_face, ...
                                      nodelist_for_face, edgelist_for_face, nodelist_for_edge);
            else
                % Use provided edge and face data
                nedge = my_nedge;
                edgelist_for_face = my_edgelist_for_face;
                nodelist_for_edge = my_nodelist_for_edge;
                norm_ea_face = my_norm_ea_face;
                facelist_for_edge = my_facelist_for_edge;
                nedge_ea_node = my_nedge_ea_node;
                edgelist_for_node = my_edgelist_for_node;
            end
        
            % Sort distances
            [ds_ordered, order] = sort(ds_ea_node);
            mynodelist_ea_ds = order;
        
            % Group nodes by distance
            mynnode_ea_ds = zeros(nnode, 1);
            ds_last = ds_ordered(1);
            mynnode_ea_ds(1) = 1;
            ndistance = 1;
        
            for i = 2:nnode
                if ds_ordered(i) > ds_last + dzmin
                    ndistance = ndistance + 1;
                    mynnode_ea_ds(ndistance) = 1;
                    ds_last = ds_ordered(i);
                else
                    mynnode_ea_ds(ndistance) = mynnode_ea_ds(ndistance) + 1;
                end
            end
        
            % Initialize lower and upper bounds
            ds_lower = ds_ordered(1);
            vf_lower = 0.0;
            ds_upper = ds_ordered(2);
        
            % Iterate through distances to find bounds
            for idx = 2:ndistance
                ds = ds_ordered(idx);
        
                % Classify nodes relative to the current z-level
                node_loc(:) = 0;
                for n = 1:nnode
                    if coords(n, 3) < ds - dzmin
                        node_loc(n) = -1; % Below
                    elseif coords(n, 3) > ds + dzmin
                        node_loc(n) = 1; % Above
                    end
                end
        
                % Find interface and divide polyhedron
                [nnode_all, coords_all, nnode_interface, nodelist_for_interface, ...
                 interface_node_to_edge] = obj.find_interface_z(nface, nedge, nnode, coords, ...
                                                            nnode_ea_face, nodelist_for_face, ...
                                                            edgelist_for_face, nodelist_for_edge, ...
                                                            nedge_ea_node, edgelist_for_node, ...
                                                            facelist_for_edge, node_loc, ds);
        
                [nnode_tot, coords_tot, nface_lower, nnode_ea_face_lower, ...
                 nodelist_for_face_lower, face2old_lower, nface_upper, ...
                 nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper, ...
                 faceindexlow_of_upper_interface, faceindexhgh_of_upper_interface] = ...
                 obj.divide_z(nface, nnode, coords, ds, node_loc, nnode_ea_face, ...
                          nodelist_for_face, nedge, edgelist_for_face, nodelist_for_edge, ...
                          nedge_ea_node, edgelist_for_node, facelist_for_edge, ...
                          norm_ea_face, nnode_all, coords_all, nnode_interface, ...
                          nodelist_for_interface, interface_node_to_edge);
        
                % Calculate volume of lower polyhedron
                if nface_lower > 0
                    vol = obj.cal_vol(nface_lower, nnode_tot, coords_tot, ...
                                  nnode_ea_face_lower, nodelist_for_face_lower);
                    vf = vol / volume;
                end
        
                % Check bounds
                if abs(vf - vf_to_match) < accuracy
                    vol_matched = 1;
                    ds_lower = ds;
                    ds_upper = ds;
                    vf_lower = vf;
                    vf_upper = vf;
                    break;
                elseif vf < vf_to_match
                    ds_lower = ds;
                    vf_lower = vf;
                else
                    ds_upper = ds;
                    vf_upper = vf;
                    break;
                end
            end
        
            % Intersection details
            [ifintersected_face, nodelist_intersected_face, nodes2_in_between] = ...
                obj.interface_face_intersect(nface, nedge, nnode, nnode_ea_face, ...
                                         nodelist_for_face, edgelist_for_face, ...
                                         nodelist_for_edge, facelist_for_edge, ...
                                         nedge_ea_node, edgelist_for_node, ...
                                         nnode_tot, nnode_interface, ...
                                         nodelist_upper_interface, interface_node_to_edge);
        
            % Cleanup temporary data
            if my_nedge <= 0
                clear norm_ea_face edgelist_for_face nodelist_for_edge facelist_for_edge;
                clear nedge_ea_node edgelist_for_node;
            end
        end

        function [distance, nnode_output, coords_output, nnode_interface, ...
                  nodelist_for_interface, faceindexlow_of_interface, ...
                  faceindexhgh_of_interface, ifintersected_face, ...
                  nodelist_intersected_face, nodes2_in_between, ...
                  nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, ...
                  face2old_lower, nface_upper, nnode_ea_face_upper, ...
                  nodelist_for_face_upper, face2old_upper] = ...
                  cal_distance_z(obj, vf_to_match, volume, nface, nnode, coords, ...
                                 nnode_ea_face, nodelist_for_face, ...
                                 node_order_for_ds, ds_ea_node)
                
            % cal_distance_z - Calculate the distance along the z-axis for a material interface.
            %
            % This function computes the distance along the z-axis (`distance`) that defines 
            % a material interface for a 3D polyhedron. It calculates the interface geometry 
            % and partitions the polyhedron into upper and lower regions.
            %
            % Inputs:
            %   vf_to_match                  - (double) Target volume fraction to match.
            %   volume                       - (double) Total volume of the polyhedron.
            %   nface                        - (int) Number of faces in the polyhedral zone.
            %   nnode                        - (int) Number of nodes in the polyhedral zone.
            %   coords                       - (double array) Coordinates of all nodes, 
            %                                  size: [nnode, 3].
            %   nnode_ea_face                - (int array) Number of nodes per face, 
            %                                  size: [nface, 1].
            %   nodelist_for_face            - (int array) Node indices for each face, 
            %                                  size: [sum(nnode_ea_face), 1].
            %   node_order_for_ds            - (int array) Order of nodes based on their 
            %                                  distance values, size: [nnode, 1].
            %   ds_ea_node                   - (double array) Distance values for each node, 
            %                                  size: [nnode, 1].
            %
            % Outputs:
            %   distance                     - (double) Calculated z-distance for the material 
            %                                  interface.
            %   nnode_output                 - (int) Total number of nodes, including interface 
            %                                  nodes.
            %   coords_output                - (double array) Coordinates of all nodes after 
            %                                  including interface nodes, size: [nnode_output, 3].
            %   nnode_interface              - (int) Number of nodes on the material interface.
            %   nodelist_for_interface       - (int array) Node indices on the material 
            %                                  interface, size: [nnode_interface, 1].
            %   faceindexlow_of_interface    - (int) Index of the interface face in the lower 
            %                                  polyhedron.
            %   faceindexhgh_of_interface    - (int) Index of the interface face in the upper 
            %                                  polyhedron.
            %   ifintersected_face           - (int array) Flags indicating whether each face 
            %                                  intersects with the interface, size: [nface, 1].
            %   nodelist_intersected_face    - (int array) Node indices where the interface 
            %                                  intersects each face, size: [2*nface, 1].
            %   nodes2_in_between            - (int array) Contains the two original nodes 
            %                                  surrounding each interface node, 
            %                                  size: [nnode_interface, 2].
            %   nface_lower                  - (int) Number of faces in the lower polyhedron.
            %   nnode_ea_face_lower          - (int array) Number of nodes per face in the 
            %                                  lower polyhedron, size: [nface_lower, 1].
            %   nodelist_for_face_lower      - (int array) Node indices for each face in the 
            %                                  lower polyhedron, size: [sum(nnode_ea_face_lower), 1].
            %   face2old_lower               - (int array) Mapping from lower polyhedron 
            %                                  faces to original mesh faces, size: [nface_lower, 1].
            %   nface_upper                  - (int) Number of faces in the upper polyhedron.
            %   nnode_ea_face_upper          - (int array) Number of nodes per face in the 
            %                                  upper polyhedron, size: [nface_upper, 1].
            %   nodelist_for_face_upper      - (int array) Node indices for each face in the 
            %                                  upper polyhedron, size: [sum(nnode_ea_face_upper), 1].
            %   face2old_upper               - (int array) Mapping from upper polyhedron 
            %                                  faces to original mesh faces, size: [nface_upper, 1].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   distance                     - Updated to the calculated z-distance.
            %   nnode_output, coords_output  - Updated with total node count and coordinates 
            %                                  after adding interface nodes.
            %   nnode_interface, nodelist_for_interface - Updated with interface node 
            %                                             information.
            %   faceindexlow_of_interface, faceindexhgh_of_interface - Indices of the interface 
            %                                                         face in the lower and 
            %                                                         upper polyhedrons.
            %   ifintersected_face, nodelist_intersected_face, nodes2_in_between - Updated with 
            %                                                                     intersection details.
            %   nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower - 
            %                                        Define the lower polyhedron.
            %   nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper - 
            %                                        Define the upper polyhedron.
            %
            % Original C declaration:
            % void cal_distance_z(double vf_to_match, double volume,
            %                     int nface, int nnode,
            %                     double *coords,
            %                     int *nnode_ea_face, int *nodelist_for_face,
            %                     int *node_order_for_ds, double *ds_ea_node,
            %                     double *distance,
            %                     int *nnode_output, double **coords_output,
            %                     int *nnode_interface, int **nodelist_for_interface,
            %                     int *faceindexlow_of_interface, int *faceindexhgh_of_interface,
            %                     int *ifintersected_face, int *nodelist_intersected_face, 
            %                     int **nodes2_in_between,
            %                     int *nface_lower, int *nnode_ea_face_lower, 
            %                     int *nodelist_for_face_lower, int *face2old_lower,
            %                     int *nface_upper, int *nnode_ea_face_upper, 
            %                     int *nodelist_for_face_upper, int *face2old_upper);
        
            % Constants and initialization
            dim = 3;
            distance = 0.0;
            nnode_interface = 0;
            nodelist_for_interface = [];
            nnode_output = 0;
            coords_output = [];
            ifintersected_face = zeros(nface, 1);
            nodelist_intersected_face = zeros(2 * nface, 1);
            nodes2_in_between = [];
            nface_lower = 0;
            nnode_ea_face_lower = [];
            nodelist_for_face_lower = [];
            face2old_lower = [];
            nface_upper = 0;
            nnode_ea_face_upper = [];
            nodelist_for_face_upper = [];
            face2old_upper = [];
        
            % Allocate memory for normals and edge data
            norm_ea_face = zeros(nface, dim);
            [nedge, edgelist_for_face, nodelist_for_edge] = ...
                obj.get_edges(nface, nnode, nnode_ea_face, nodelist_for_face);
        
            facelist_for_edge = zeros(2 * nedge, 1);
            [nedge_ea_node, edgelist_for_node] = ...
                obj.get_faces_ea_edge(nface, nedge, nnode, nnode_ea_face, ...
                                  nodelist_for_face, edgelist_for_face, ...
                                  nodelist_for_edge, facelist_for_edge);
        
            norm_ea_face = obj.cal_norm_of_face(nface, nnode, coords, ...
                                            nnode_ea_face, nodelist_for_face, ...
                                            edgelist_for_face);
        
            % Compute bounds for distance
            [ds_lower, ds_upper, vf_lower, vf_upper, nnode_output, coords_output, ...
             if_matched, nnode_upper_interface, nodelist_upper_interface, ...
             faceindexlow_of_interface, faceindexhgh_of_interface, ...
             ifintersected_face, nodelist_intersected_face, nodes2_in_between, ...
             nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, ...
             face2old_lower, nface_upper, nnode_ea_face_upper, ...
             nodelist_for_face_upper, face2old_upper] = ...
                obj.bounds_z(vf_to_match, volume, nface, nnode, coords, ...
                         nnode_ea_face, nodelist_for_face, node_order_for_ds, ...
                         ds_ea_node, nedge, edgelist_for_face, ...
                         nodelist_for_edge, nedge_ea_node, edgelist_for_node, ...
                         facelist_for_edge, norm_ea_face);
        
            if if_matched
                % Bounds calculation succeeded
                distance = ds_upper;
                nnode_interface = nnode_upper_interface;
                nodelist_for_interface = nodelist_upper_interface;
            else
                % Perform bisection to refine distance
                [distance, nnode_output, coords_output, nnode_interface, ...
                 nodelist_for_interface, interface_node_to_edge, ...
                 faceindexlow_of_interface, faceindexhgh_of_interface, ...
                 nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, ...
                 face2old_lower, nface_upper, nnode_ea_face_upper, ...
                 nodelist_for_face_upper, face2old_upper] = ...
                    obj.bisection_z(vf_to_match, volume, ds_lower, ds_upper, vf_lower, ...
                                vf_upper, nface, nedge, nnode, coords, ...
                                nnode_ea_face, nodelist_for_face, ...
                                edgelist_for_face, nodelist_for_edge);
            end
        
            % Clean up temporary resources
            clear norm_ea_face facelist_for_edge edgelist_for_face nodelist_for_edge;
            clear nedge_ea_node edgelist_for_node;
        end
 
        function [distance, nnode_output, coords_output, nnode_interface, ...
                  nodelist_for_interface, faceindexlow_of_interface, ...
                  faceindexhgh_of_interface, ifintersected_face, ...
                  nodelist_intersected_face, nodes2_in_between, ...
                  nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, ...
                  face2old_lower, nface_upper, nnode_ea_face_upper, ...
                  nodelist_for_face_upper, face2old_upper] = ...
                  cal_distance(obj, vf_to_match, volume, norm, nface, nnode, ...
                               coords, nnode_ea_face, nodelist_for_face)
            % cal_distance - Calculate the distance of a material interface in 3D.
            %
            % This function computes the distance (`distance`) of a material interface 
            % in 3D space, defined by a plane with the equation:
            % `norm[0] * x + norm[1] * y + norm[2] * z = distance`.
            % It also partitions the polyhedron into upper and lower regions, calculates 
            % the interface geometry, and provides additional intersection details.
            %
            % Inputs:
            %   vf_to_match                  - (double) Target volume fraction to match.
            %   volume                       - (double) Total volume of the polyhedron.
            %   norm                         - (double array) Normal vector defining the 
            %                                  material interface, size: [3, 1].
            %   nface                        - (int) Number of faces in the polyhedral zone.
            %   nnode                        - (int) Number of nodes in the polyhedral zone.
            %   coords                       - (double array) Coordinates of all nodes, 
            %                                  size: [nnode, 3].
            %   nnode_ea_face                - (int array) Number of nodes per face, 
            %                                  size: [nface, 1].
            %   nodelist_for_face            - (int array) Node indices for each face, 
            %                                  size: [sum(nnode_ea_face), 1].
            %
            % Outputs:
            %   distance                     - (double) Calculated distance for the material 
            %                                  interface.
            %   nnode_output                 - (int) Total number of nodes, including interface 
            %                                  nodes.
            %   coords_output                - (double array) Coordinates of all nodes after 
            %                                  adding interface nodes, size: [nnode_output, 3].
            %   nnode_interface              - (int) Number of nodes on the material interface.
            %   nodelist_for_interface       - (int array) Node indices on the material 
            %                                  interface, size: [nnode_interface, 1].
            %   faceindexlow_of_interface    - (int) Index of the interface face in the lower 
            %                                  polyhedron.
            %   faceindexhgh_of_interface    - (int) Index of the interface face in the upper 
            %                                  polyhedron.
            %   ifintersected_face           - (int array) Flags indicating whether each face 
            %                                  intersects with the interface, size: [nface, 1].
            %   nodelist_intersected_face    - (int array) Node indices where the interface 
            %                                  intersects each face, size: [2*nface, 1].
            %   nodes2_in_between            - (int array) Contains the two original nodes 
            %                                  surrounding each interface node, 
            %                                  size: [nnode_interface, 2].
            %   nface_lower                  - (int) Number of faces in the lower polyhedron.
            %   nnode_ea_face_lower          - (int array) Number of nodes per face in the 
            %                                  lower polyhedron, size: [nface_lower, 1].
            %   nodelist_for_face_lower      - (int array) Node indices for each face in the 
            %                                  lower polyhedron, size: [sum(nnode_ea_face_lower), 1].
            %   face2old_lower               - (int array) Mapping from lower polyhedron 
            %                                  faces to original mesh faces, size: [nface_lower, 1].
            %   nface_upper                  - (int) Number of faces in the upper polyhedron.
            %   nnode_ea_face_upper          - (int array) Number of nodes per face in the 
            %                                  upper polyhedron, size: [nface_upper, 1].
            %   nodelist_for_face_upper      - (int array) Node indices for each face in the 
            %                                  upper polyhedron, size: [sum(nnode_ea_face_upper), 1].
            %   face2old_upper               - (int array) Mapping from upper polyhedron 
            %                                  faces to original mesh faces, size: [nface_upper, 1].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   distance                     - Updated to the calculated distance.
            %   nnode_output, coords_output  - Updated with total node count and coordinates 
            %                                  after adding interface nodes.
            %   nnode_interface, nodelist_for_interface - Updated with interface node 
            %                                             information.
            %   faceindexlow_of_interface, faceindexhgh_of_interface - Indices of the interface 
            %                                                         face in the lower and 
            %                                                         upper polyhedrons.
            %   ifintersected_face, nodelist_intersected_face, nodes2_in_between - Updated with 
            %                                                                     intersection details.
            %   nface_lower, nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower - 
            %                                        Define the lower polyhedron.
            %   nface_upper, nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper - 
            %                                        Define the upper polyhedron.
            %
            % Original C declaration:
            % void cal_distance(double vf_to_match, double volume, double *norm,
            %                   int nface, int nnode, double *coords,
            %                   int *nnode_ea_face, int *nodelist_for_face,
            %                   double *distance,
            %                   int *nnode_output, double **coords_output,
            %                   int *nnode_interface, int **nodelist_for_interface,
            %                   int *faceindexlow_of_interface, int *faceindexhgh_of_interface,
            %                   int *ifintersected_face, int *nodelist_intersected_face, 
            %                   int **nodes2_in_between,
            %                   int *nface_lower, int *nnode_ea_face_lower, 
            %                   int *nodelist_for_face_lower, int *face2old_lower,
            %                   int *nface_upper, int *nnode_ea_face_upper, 
            %                   int *nodelist_for_face_upper, int *face2old_upper);
        
            % Initialize variables
            dim = 3;
            distance = 0.0;
            coords_rotated = zeros(nnode, dim);
            coords_output = [];
            nnode_output = 0;
            zmin = inf;
            zmax = -inf;
        
            % Calculate zmin and zmax
            for i = 1:nnode
                zmin = min(zmin, coords(i, 3));
                zmax = max(zmax, coords(i, 3));
            end
            dzmin = 1.0e-6 * (zmax - zmin);
        
            % Rotate coordinates to align with the interface normal
            coords_rotated = obj.rotate_to_norm(norm, nnode, coords);
        
            % Set the rotated normal
            norm_rotated = [0.0, 0.0, 1.0];
        
            % Initialize arrays for sorting nodes
            node_order_for_ds = zeros(nnode, 1);
            ds_ea_node = zeros(nnode, 1);
        
            % Order nodes along the rotated normal
            [node_order_for_ds, ds_ea_node] = obj.order_nodes_along_norm(dim, norm_rotated, ...
                                                                     nnode, coords_rotated);
        
            % Calculate distance along the z-axis
            [distance, nnode_tot, coords_tot, nnode_interface, nodelist_for_interface, ...
             faceindexlow_of_interface, faceindexhgh_of_interface, ifintersected_face, ...
             nodelist_intersected_face, nodes2_in_between, nface_lower, ...
             nnode_ea_face_lower, nodelist_for_face_lower, face2old_lower, nface_upper, ...
             nnode_ea_face_upper, nodelist_for_face_upper, face2old_upper] = ...
                obj.cal_distance_z(vf_to_match, volume, nface, nnode, coords_rotated, ...
                               nnode_ea_face, nodelist_for_face, node_order_for_ds, ds_ea_node);
        
            % Rotate the output coordinates back to the original frame
            coords_output = obj.rotate_back(norm, nnode_tot, coords_tot);
        
            % Assign outputs
            nnode_output = nnode_tot;
        
            % Cleanup
            clear coords_tot coords_rotated node_order_for_ds ds_ea_node;
        end
        
        function [ifinterface, nnode_new, coords_new, nnode_for_interface, ...
                  nodelist_for_interface, faceindexlow_of_interface, ...
                  faceindexhgh_of_interface, ifintersected_face, ...
                  nodelist_intersected_face, nodes2_in_between, ...
                  nface_lower, nnode_for_face_lower, nodelist_for_face_lower, ...
                  face2old_lower, nface_upper, nnode_for_face_upper, ...
                  nodelist_for_face_upper, face2old_upper] = ...
                  interface3d(obj, x0_scale, dx_scale, nface, nnode, ...
                              coords, nnode_for_face, nodelist_for_face, ...
                              cell_vol, vf_to_match, normal)
            % interface3d - Compute the material interface for a 3D polyhedron.
            %
            % This function calculates the material interface for a 3D polyhedron based 
            % on the provided volume fraction and normal vector. The interface is represented 
            % as a plane that divides the polyhedron into lower and upper regions. It also 
            % computes additional intersection details.
            %
            % Inputs:
            %   x0_scale                  - (double array) Scaling origin, size: [3, 1].
            %   dx_scale                  - (double) Scaling factor.
            %   nface                     - (int) Number of faces in the polyhedron.
            %   nnode                     - (int) Number of nodes in the polyhedron.
            %   coords                    - (double array) Coordinates of all nodes, 
            %                               size: [nnode, 3].
            %   nnode_for_face            - (int array) Number of nodes per face, 
            %                               size: [nface, 1].
            %   nodelist_for_face         - (int array) Node indices for each face, 
            %                               size: [sum(nnode_for_face), 1].
            %   cell_vol                  - (double) Total volume of the polyhedron.
            %   vf_to_match               - (double) Target volume fraction to match.
            %   normal                    - (double array) Normal vector defining the 
            %                               material interface, size: [3, 1].
            %
            % Outputs:
            %   ifinterface               - (int) Flag indicating if an interface was found 
            %                               (1 if found, 0 otherwise).
            %   nnode_new                 - (int) Total number of nodes, including interface 
            %                               nodes.
            %   coords_new                - (double array) Coordinates of all nodes after 
            %                               including interface nodes, size: [nnode_new, 3].
            %   nnode_for_interface       - (int) Number of nodes on the material interface.
            %   nodelist_for_interface    - (int array) Node indices on the material 
            %                               interface, size: [nnode_for_interface, 1].
            %   faceindexlow_of_interface - (int) Index of the interface face in the lower 
            %                               polyhedron.
            %   faceindexhgh_of_interface - (int) Index of the interface face in the upper 
            %                               polyhedron.
            %   ifintersected_face        - (int array) Flags indicating whether each face 
            %                               intersects with the interface, size: [nface, 1].
            %   nodelist_intersected_face - (int array) Node indices where the interface 
            %                               intersects each face, size: [2*nface, 1].
            %   nodes2_in_between         - (int array) Contains the two original nodes 
            %                               surrounding each interface node, 
            %                               size: [nnode_for_interface, 2].
            %   nface_lower               - (int) Number of faces in the lower polyhedron.
            %   nnode_for_face_lower      - (int array) Number of nodes per face in the 
            %                               lower polyhedron, size: [nface_lower, 1].
            %   nodelist_for_face_lower   - (int array) Node indices for each face in the 
            %                               lower polyhedron, size: [sum(nnode_for_face_lower), 1].
            %   face2old_lower            - (int array) Mapping from lower polyhedron 
            %                               faces to original mesh faces, size: [nface_lower, 1].
            %   nface_upper               - (int) Number of faces in the upper polyhedron.
            %   nnode_for_face_upper      - (int array) Number of nodes per face in the 
            %                               upper polyhedron, size: [nface_upper, 1].
            %   nodelist_for_face_upper   - (int array) Node indices for each face in the 
            %                               upper polyhedron, size: [sum(nnode_for_face_upper), 1].
            %   face2old_upper            - (int array) Mapping from upper polyhedron 
            %                               faces to original mesh faces, size: [nface_upper, 1].
            %
            % Uses:
            %   None.
            %
            % Modifies:
            %   ifinterface               - Indicates if the interface was found.
            %   nnode_new, coords_new     - Updated with total node count and coordinates 
            %                               after including interface nodes.
            %   nnode_for_interface, nodelist_for_interface - Updated with interface node 
            %                                                 information.
            %   faceindexlow_of_interface, faceindexhgh_of_interface - Indices of the interface 
            %                                                         face in the lower and 
            %                                                         upper polyhedrons.
            %   ifintersected_face, nodelist_intersected_face, nodes2_in_between - Updated with 
            %                                                                     intersection details.
            %   nface_lower, nnode_for_face_lower, nodelist_for_face_lower, face2old_lower - 
            %                                        Define the lower polyhedron.
            %   nface_upper, nnode_for_face_upper, nodelist_for_face_upper, face2old_upper - 
            %                                        Define the upper polyhedron.
            %
            % Original C declaration:
            % void interface3d(double *x0_scale, double dx_scale,
            %                  int nface, int nnode, double *coords,
            %                  int *nnode_for_face, int *nodelist_for_face,
            %                  double cell_vol, double vf_to_match, double *normal,
            %                  int *ifinterface,
            %                  int *nnode_new,  double **coords_new,
            %                  int *nnode_for_interface, int **nodelist_for_interface,
            %                  int *faceindexlow_of_interface, int *faceindexhgh_of_interface,
            %                  int *ifintersected_face, int *nodelist_intersected_face, 
            %                  int **nodes2_in_between,
            %                  int *nface_lower, int *nnode_for_face_lower, 
            %                  int *nodelist_for_face_lower, int *face2old_lower,
            %                  int *nface_upper, int *nnode_for_face_upper, 
            %                  int *nodelist_for_face_upper, int *face2old_upper);
        
        
            % Initialize variables
            dim = 3;
            coords_scaled = [];
            volume = cell_vol;
            distance = 0.0;
        
            % Scale coordinates if scaling parameters are provided
            if ~isempty(x0_scale)
                coords_scaled = obj.coords_scale(dim, x0_scale, dx_scale, nnode, coords);
                volume = cell_vol / (dx_scale^3);
            else
                coords_scaled = coords;
            end
        
            % Calculate the interface distance and other properties
            [distance, nnode_new, coords_new, nnode_for_interface, nodelist_for_interface, ...
             faceindexlow_of_interface, faceindexhgh_of_interface, ifintersected_face, ...
             nodelist_intersected_face, nodes2_in_between, nface_lower, ...
             nnode_for_face_lower, nodelist_for_face_lower, face2old_lower, nface_upper, ...
             nnode_for_face_upper, nodelist_for_face_upper, face2old_upper] = ...
                obj.cal_distance(vf_to_match, volume, normal, nface, nnode, coords_scaled, ...
                                 nnode_for_face, nodelist_for_face);
        
            % Scale back the coordinates if scaling was applied
            if ~isempty(x0_scale)
                coords_new = obj.coords_scale_back(dim, x0_scale, dx_scale, nnode_new, coords_new);
            end
        
            % Indicate if the interface was found
            ifinterface = ~isempty(nnode_for_interface) && nnode_for_interface > 0;
        end

        function [nnode_final, coords_final, nnode_for_minterface, nodes_for_minterface, ...
          nface_for_mat, nnode_for_face_ea_mat, nodelist_for_face_ea_mat] = ...
          reconstruct3d_nmat_pagosa(~, xl, dx, nmat_mesh, matid_mesh, vf_mesh)
            % reconstruct3d_nmat_pagosa - Reconstruct the 3D geometry for multiple materials
            %                              in a PAGOSA-style mesh.
            %
            % This function reconstructs the 3D geometry of a cell that contains multiple
            % materials using the PAGOSA method. It calculates the material interfaces and
            % reconstructs the polyhedra representing each material within the cell.
            %
            % Inputs:
            %   xl                  - (double array) The lower bounds of the cell in the 
            %                         x, y, and z directions, size: [3, 1].
            %   dx                  - (double array) The grid spacing in the x, y, and z 
            %                         directions, size: [3, 1].
            %   nmat_mesh           - (int array) The number of materials in each cell of 
            %                         the mesh, size: [nx, ny, nz].
            %   matid_mesh          - (cell array of int arrays) The material IDs in each 
            %                         cell of the mesh, with each cell containing an int 
            %                         array of size [nmat_cell, 1].
            %   vf_mesh             - (cell array of double arrays) The volume fractions of 
            %                         each material in each cell of the mesh, with each 
            %                         cell containing a double array of size [nmat_cell, 1].
            %
            % Outputs:
            %   nnode_final         - (int) The total number of nodes in the final reconstructed
            %                         geometry.
            %   coords_final        - (double array) The final coordinates of the nodes after 
            %                         reconstruction, size: [nnode_final, 3].
            %   nnode_for_minterface - (int array) The number of nodes for each material 
            %                         interface, size: [nmat, 1].
            %   nodes_for_minterface - (cell array of int arrays) The list of node indices for 
            %                         each material interface, with each element of the cell array
            %                         containing an int array.
            %   nface_for_mat       - (int array) The number of faces for each material 
            %                         polyhedron, size: [nmat, 1].
            %   nnode_for_face_ea_mat - (cell array of int arrays) The number of nodes for each 
            %                         face of each material polyhedron, with each element of the 
            %                         cell array containing an int array.
            %   nodelist_for_face_ea_mat - (cell array of int arrays) The list of node indices 
            %                         for each face of each material polyhedron, with each 
            %                         element of the cell array containing an int array.
            %
            % Original C declaration:
            % void reconstruct3d_nmat_pagosa(double *xl, double *dx,
            %                                int ***nmat_mesh, int ****matid_mesh, double ****vf_mesh,
            %                                int *nnode_final, double **coords_final,
            %                                int *nnode_for_minterface, int ***nodes_for_minterface,
            %                                int *nface_for_mat, int ***nnode_for_face_ea_mat,
            %                                int ***nodelist_for_face_ea_mat);
        
            global mesh tiny small vof3d_lib;

            % Initialize constants and parameters
            dim = 3;
            nn_ea_face = 4;
            dx_max = max(dx);
            mydx = dx / dx_max;
            cell_vol = prod(mydx);
        
            % Initialize coordinates
            nnode = 8;
            coords = zeros(nnode, dim);
            coords(1, :) = [0, 0, 0];
            coords(2, :) = [mydx(1), 0, 0];
            coords(3, :) = [mydx(1), mydx(2), 0];
            coords(4, :) = [0, mydx(2), 0];
            coords(5, :) = [0, 0, mydx(3)];
            coords(6, :) = [mydx(1), 0, mydx(3)];
            coords(7, :) = [mydx(1), mydx(2), mydx(3)];
            coords(8, :) = [0, mydx(2), mydx(3)];
        
            % Initialize material and node information
            nmat = nmat_mesh(2, 2, 2);  % Extract the number of materials in the central cell
            nnode_max = nnode * nmat;
            coords_final = zeros(nnode_max, dim);
            coords_final(1:nnode, :) = coords;
            nnode_final = nnode;
            
            % Initialize face information
            nface = 6;
            nnode_for_face = repmat(4, nface, 1);
            
            % Face 0 (x0)
            nodelist_for_face = zeros(nface * nn_ea_face, 1);
            nodelist_for_face(1:4) = [1, 5, 8, 4];
            
            % Face 1 (x1)
            offset = nn_ea_face;
            nodelist_for_face(offset + (1:4)) = [2, 3, 7, 6];
            
            % Face 2 (y0)
            offset = 2 * nn_ea_face;
            nodelist_for_face(offset + (1:4)) = [1, 2, 6, 5];
            
            % Face 3 (y1)
            offset = 3 * nn_ea_face;
            nodelist_for_face(offset + (1:4)) = [3, 4, 8, 7];
            
            % Face 4 (z0)
            offset = 4 * nn_ea_face;
            nodelist_for_face(offset + (1:4)) = [1, 4, 3, 2];
            
            % Face 5 (z1)
            offset = 5 * nn_ea_face;
            nodelist_for_face(offset + (1:4)) = [5, 6, 7, 8];
            
            % Allocate arrays for intersected faces and material face indices
            ifintersected_face = zeros(nface, 1);
            nodelist_intersected_face = zeros(nface * 2, 1);
            face2old_lower = zeros(nface + 1, 1);
            face2old_upper = zeros(nface + 1, 1);

            % Allocate output
            nface_for_mat = zeros(nmat,1);
        
            % Initialize temporary storage for material interface and face data
            nodes_for_minterface_tmp = cell(nmat-1,1);   % zeros(nmat - 1, 8);
            nnode_for_face_ea_mat_tmp = cell(nmat,1);    % zeros(nmat, nface + 1);
            nodelist_for_face_ea_mat_tmp = cell(nmat,1); % zeros(nmat, (nface + 1)*8);
            
            % Initialize volume fractions (vfs)
            vfs = zeros(2, 3, 3, 3);
            
            % Update volume fractions for materials
            for m = 1:(nmat - 1)
                for k = 1:3
                    for j = 1:3
                        for i = 1:3
                            if k == 2 && j == 2 && i == 2
                                continue;
                            end
                            vfs(2, k, j, i) = 0.0; % Outer material fractions initialized
                        end
                    end
                end
                vfs(1, 2, 2, 2) = vfs(1, 2, 2, 2) + vf_mesh(2, 2, 2, m); % Inner material
                vfs(2, 2, 2, 2) = 1.0 - vfs(1, 2, 2, 2); % Outer material fraction
            
                ids = matid_mesh(2, 2, 2, 1:nmat);
                for k = 1:3
                    for j = 1:3
                        for i = 1:3
                            if i == 2 && j == 2 && k == 2
                                continue;
                            end
            
                            nm_neighb = nmat_mesh(k, j, i);
                            ids_neighb = matid_mesh(k, j, i, 1:nmat);
                            for m_neighb = 1:nm_neighb
                                if ids_neighb(m_neighb) == ids(m)
                                    vfs(1, k, j, i) = vfs(1, k, j, i) + vf_mesh(k, j, i, m_neighb);
                                else
                                    for m_other = (m + 1):nmat
                                        if ids_neighb(m_neighb) == ids(m_other)
                                            vfs(2, k, j, i) = vfs(2, k, j, i) + vf_mesh(k, j, i, m_neighb);
                                            break;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
        
                % Calculate gradients and determine the interface normal
                priority_max = 0.0;
                idx = 0;
                norms = zeros(2, dim); % Normals for materials
                grad = zeros(1, dim);  % Gradient vector
                
                for m_other = 1:2
                    mynorm = norms(m_other, :);
                    grad = mesh.cal_cell_zgrad3d(dx, squeeze(vfs(m_other, :, :, :)));
                    vsq = sum(grad .* grad);
                    mynorm = -grad; % Negate the gradient for normal calculation
                    priority = sqrt(vsq) * sqrt(vfs(m_other, 2, 2, 2));
                    
                    if priority > priority_max
                        idx = m_other;
                        priority_max = priority;
                    end
                    
                    % Normalize the normal
                    factor = 1.0 / (sqrt(vsq) + tiny);
                    norms(m_other, :) = mynorm * factor;
                end
                
                % Select the final normal vector
                mynorm = norms(idx, :);
                if idx == 2
                    norm = -mynorm;
                else
                    norm = mynorm;
                end

                %%%% TODO: add "check norm"
                factor = sum(norm.^2);
                if factor < small
                    norm(1) = 1.0;
                end
                
                % Initialize variables for interface reconstruction
                nnode_out = 0;
                coords_out = [];
                nnode_interface = 0;
                nodelist_interface = [];
                
                nface_lower = 0;
                nnode_for_face_lower = zeros(nface + 1, 1);
                nodelist_for_face_lower = zeros((nface + 1) * 8, 1);
                
                nface_upper = 0;
                nnode_for_face_upper = zeros(nface + 1, 1);
                nodelist_for_face_upper = zeros((nface + 1) * 8, 1);
                
                nodes2_in_between = [];
                
                % Perform 3D interface reconstruction for the current material
                % ===============> use C call!
                %[ifinterface, nnode_out, coords_out, nnode_interface, nodelist_interface, ...
                % faceindexlow, faceindexhgh, ifintersected_face, nodelist_intersected_face, ...
                % nodes2_in_between, nface_lower, nnode_for_face_lower, ...
                % nodelist_for_face_lower, face2old_lower, nface_upper, ...
                % nnode_for_face_upper, nodelist_for_face_upper, face2old_upper] = ...
                % obj.interface3d([], dx_max, nface, nnode, coords, ...
                % nnode_for_face, nodelist_for_face, cell_vol, vfs(1, 2, 2, 2), ...
                % norm);
                [coords_out, nodelist_interface, ~, ~,...
                nface_lower, nnode_for_face_lower, nodelist_for_face_lower, ... 
                nface_upper, nnode_for_face_upper, nodelist_for_face_upper] = ...
                vof3d_lib.interface3d(dx_max, nface, nnode, coords, ...
                            nnode_for_face, nodelist_for_face, cell_vol, ...
                            vfs(1,2,2,2), norm);
                nnode_out = size(coords_out, 1)/3;
                coords_out = reshape(coords_out, 3, nnode_out)';
                
                % Check and resize coordinates array if needed
                if nnode_final + (nnode_out - nnode) > nnode_max
                    nnode_max = nnode_max + nmat * (nnode_out - nnode);
                    coords_final = [coords_final; zeros(nnode_max * dim - size(coords_final, 1), dim)];
                end
                
                % Add new coordinates to the final coordinates array
                nn_new = nnode_out - nnode;
                coords_final(nnode_final+1:nnode_final+nn_new, :) = coords_out(nnode+1:nnode_out,:);
                
                % Free memory for temporary output coordinates
                clear coords_out;
                
                % Store the number of nodes on the material interface
                nnode_for_minterface(m) = nnode_interface;
                
                % Adjust node indices for the interface nodes
                if nnode_interface > 2
                    for i = 1:nnode_interface
                        n = nodelist_interface(i);
                        if n >= nnode
                            nnew = n - nnode;
                            nodelist_interface(i) = nnode_final + nnew;
                        end
                    end
                end
                
                % Adjust node indices for the lower face nodes
                nodelist = nodelist_for_face_lower;
                for f = 1:nface_lower
                    nn = nnode_for_face_lower(f);
                    for i = 1:nn
                        n = nodelist(i);
                        if n > nnode
                            nnew = n - nnode;
                            nodelist(i) = nnode_final + nnew;
                        end
                    end
                    nodelist = nodelist(nn+1:end); % Move to the next face
                end
                
                % Store results for the current material
                nface_for_mat(m) = nface_lower;
                nodes_for_minterface_tmp{m} = nodelist_interface;
                nnode_for_face_ea_mat_tmp{m} = nnode_for_face_lower;
                nodelist_for_face_ea_mat_tmp{m} = nodelist_for_face_lower;
                
                % Free resources for the upper face if not the last material
                if m < nmat - 1
                    if exist('nnode_for_face_upper', 'var'), clear nnode_for_face_upper; end
                    if exist('nodelist_for_face_upper', 'var'), clear nodelist_for_face_upper; end
                end
                
                % Update the number of final nodes
                nnode_previous = nnode_final; % For the last nodelist_upper
                nnode_final = nnode_final + nn_new;
            end
            % the C trick with the loop where the variable equals end value
            % of the loop upon exit doesn't work in Matlab
            m = nmat;
            
            % Adjust node indices for the upper face nodes
            ofs = 0;
            for f = 1:nface_upper
                nn = nnode_for_face_upper(f);
                for i = 1:nn
                    n = nodelist_for_face_upper(i+ofs);
                    if n > nnode
                        nnew = n - nnode;
                        nodelist_for_face_upper(i+ofs) = nnode_previous + nnew;
                    end
                end
                ofs = ofs + nn;
            end
            
            nface_for_mat(m) = nface_upper;
            nnode_for_face_ea_mat_tmp{m} = nnode_for_face_upper;
            nodelist_for_face_ea_mat_tmp{m} = nodelist_for_face_upper;
            
            % Consolidate nodes for material interfaces
            nodes_for_minterface = cell(nmat-1, 1);
            for m = 1:nmat-1
                nn = nnode_for_minterface(m);
                nodes_for_minterface{m} = nodes_for_minterface_tmp{m}(1:nn);
            end
            
            % Consolidate number of faces and nodes for material interfaces
            nfsum = sum(nface_for_mat);
            nnode_for_face_ea_mat = cell(nmat, 1);
            nodelist_for_face_ea_mat = cell(nmat, 1);
            
            % Copy face nodes for material interfaces
            for m = 1:nmat
                nf = nface_for_mat(m);
                nnode_for_face_ea_mat{m} = nnode_for_face_ea_mat_tmp{m}(1:nf);
                lsize = sum(nnode_for_face_ea_mat{m});
                nodelist_for_face_ea_mat{m} = nodelist_for_face_ea_mat_tmp{m}(1:lsize);
            end
            
            % Scale coordinates back
            for n = 1:nnode_final
                c = coords_final(n, :);
                coords_final(n, :) = xl + dx_max * c;
            end
            
            % Free temporary arrays
            clear ifintersected_face nodelist_intersected_face face2old_lower face2old_upper coords;
        
        end

    end
end
