classdef c_geom < handle
    
    properties
        % Flags
        if_decompose_concave = 0;
        if_write_debug_viz = 0;

        % Reserved for OSO file
        nfile_oso = 0;
        reg_oso = zeros(1, 100);
        nxy_oso = zeros(1, 100);
        xys_oso cell = cell(1, 100);
        fp_oso cell = cell(1, 100);

        npart_oso = zeros(1, 100);
        nnode_part_oso cell = cell(1, 100);
        factor_part_oso cell = cell(1, 100);
        coords_part_oso cell = cell(1, 100);

        nnode_grp_oso cell = cell(1, 100);
        nedge_grp_oso cell = cell(1, 100);
        nface_grp_oso cell = cell(1, 100);
        ctr_grp_oso cell = cell(1, 100);
        coords_grp_oso cell = cell(1, 100);

        nodelist_for_edge_oso cell = cell(1, 100);
        edgelist_for_face_oso cell = cell(1, 100);
        nodelist_for_face_oso cell = cell(1, 100);
        nnode_for_face_oso cell = cell(1, 100);
        normal_of_face_oso cell = cell(1, 100);

        % Reserved for oblique ellipse
        aligned_sav = 1;
        major_axis_sav = -1;
        sinphi_sav = 0.0;
        cosphi_sav = 1.0;
        sintheta_sav = 0.0;
        costheta_sav = 1.0;

        % Reserved for conic
        sinphi_conic = 0.0;
        cosphi_conic = 1.0;
        sintheta_conic = 0.0;
        costheta_conic = 1.0;
    end

    methods

%%%%%%%%%%%%%%% implemented but not tested
        
        function [mixed, vol] = rec_rec(~, ifinquiry, geop, vcell, dim, xxl, xxr, xl, dx)
            % Function: rec_rec
            % Inputs:
            %   ifinquiry - Inquiry flag to check volume fraction or mixed/clean status
            %   geop - Type of geometry (1-Cartesian, 2-cylindrical, 3-spherical)
            %   vcell - Volume of the cell
            %   dim - Dimensionality of the problem
            %   xxl - Lower coordinates of the rectangular region
            %   xxr - Upper coordinates of the rectangular region
            %   xl - Lower coordinates of the cell
            %   dx - Cell dimensions
            %
            % Outputs:
            %   mixed - Indicates if the cell intersects the rectangular region (0 = no, 1 = yes)
            %   vol - Volume of the intersection

            % Initialize variables
            small = 1.0e-10 * dx(1);
            
            xr = xl + dx;
        
            % Check if the rectangle is completely outside the cell
            outside = false;
            for k = 1:dim
                if (xxl(k) + small > xr(k)) || (xxr(k) < xl(k) + small)
                    outside = true;
                    break;
                end
            end
        
            if outside
                vol = 0.0;
                mixed = 0;
                return;
            end
        
            % Check if the rectangle is completely inside the cell
            inside = true;
            for k = 1:dim
                if ~((xl(k) + small > xxl(k)) && (xr(k) < xxr(k) + small))
                    inside = false;
                    break;
                end
            end
        
            if inside
                vol = vcell;
                mixed = 0;
                return;
            end
        
            % Inquiry only for volume fraction, not detailed intersection
            if ifinquiry == 1
                vol = small;
                mixed = 1;
                return;
            end
        
            % The cell and rectangle are mixed
            mixed = 1;
        
            if (dim == 1) || (geop == 1)
                vol = 1.0;
                for k = 1:dim
                    z0 = max(xl(k), xxl(k));
                    z1 = min(xr(k), xxr(k));
                    vol = vol * (z1 - z0);
                end
            elseif (dim == 2) && (geop == 2)
                z0 = max(xl(2), xxl(2));
                z1 = min(xr(2), xxr(2));
                h = max(0.0, z1 - z0);
                r0 = max(xl(1), xxl(1));
                r1 = min(xr(1), xxr(1));
                dr = max(0.0, r1 - r0);
                vol = h * (r1 + r0) * dr;
            elseif (dim == 1) && (geop == 2)
                vol = 1.0;
                for k = 1:dim
                    z0 = max(xl(k), xxl(k));
                    z1 = min(xr(k), xxr(k));
                    vol = vol * (z1^2 - z0^2);
                end
            elseif (dim == 1) || (geop == 3)
                vol = (4.0 / 3.0);
                for k = 1:dim
                    z0 = max(xl(k), xxl(k));
                    z1 = min(xr(k), xxr(k));
                    vol = vol * (z1^3 - z0^3);
                end
            elseif dim == 2
                if (xl(1) < 0.0) && (xr(1) > 0.0)
                    nc = 2;
                    myxl(1, 1) = 0.0; myxr(1, 1) = -xl(1);
                    myxl(2, 1) = 0.0; myxr(2, 1) = xr(1);
                    myxxl(1, 1) = 0.0; myxxr(1, 1) = min(0.0, -xxl(1));
                    myxxl(2, 1) = max(0.0, xxl(1)); myxxr(2, 1) = xxr(1);
                elseif xl(1) < 0.0
                    nc = 1;
                    myxl(1, 1) = -xr(1); myxr(1, 1) = -xl(1);
                    myxxl(1, 1) = -xxr(1); myxxr(1, 1) = -xxl(1);
                else
                    nc = 1;
                    myxl(1, 1) = xl(1); myxr(1, 1) = xr(1);
                    myxxl(1, 1) = xxl(1); myxxr(1, 1) = xxr(1);
                end
        
                for k = 1:nc
                    myxl(k, 2) = xl(2);
                    myxr(k, 2) = xr(2);
                    myxxl(k, 2) = xxl(2);
                    myxxr(k, 2) = xxr(2);
                end
        
                vol = 0.0;
                for k = 1:nc
                    cl = myxl(k, :);
                    cr = myxr(k, :);
                    ccl = myxxl(k, :);
                    ccr = myxxr(k, :);
                    r0 = max(cl(1), ccl(1));
                    r1 = min(cr(1), ccr(1));
                    z0 = max(cl(2), ccl(2));
                    z1 = min(cr(2), ccr(2));
                    dvol = (r1 + r0) * (r1 - r0) * (z1 - z0);
                    vol = vol + dvol;
                end
            end
        end
        
%%%%%%%%%%%%%%% empty stubs
        
        function [mixed, vol] = gsph_rec(obj, ifinquiry, geop, dim, ctr, r, xl, dx)
            % Function: gsph_rec
            % Inputs:
            %   ifinquiry - Inquiry flag to check volume fraction or mixed/clean status
            %   geop - Type of geometry (1-Cartesian, 2-cylindrical, 3-spherical)
            %   dim - Dimensionality of the problem
            %   ctr - Center of the sphere
            %   r - Radius of the sphere
            %   xl - Lower coordinates of the cell
            %   dx - Cell dimensions
            %
            % Outputs:
            %   mixed - Indicates if the cell intersects the sphere (0 = no, 1 = yes)
            %   vol - Volume of the intersection

            % TODO
            mixed = 0;
            vol = 1.0;

        end
        
        function [mixed, vol] = poly2d_rec(obj, ifinquiry, geop, vcell, dim, nnode, xys, x2min, x2max, xl, dx)
            % Function: poly2d_rec
            % Inputs:
            %   ifinquiry - Inquiry flag to check volume fraction or mixed/clean status
            %   geop - Type of geometry (1-Cartesian, 2-cylindrical, 3-spherical)
            %   vcell - Volume of the cell
            %   dim - Dimensionality of the problem
            %   nnode - Number of nodes defining the polygon
            %   xys - Coordinates of the polygon vertices
            %   x2min - Minimum x2 coordinate
            %   x2max - Maximum x2 coordinate
            %   xl - Lower coordinates of the cell
            %   dx - Cell dimensions
            %
            % Outputs:
            %   mixed - Indicates if the cell intersects the polygon (0 = no, 1 = yes)
            %   vol - Volume of the intersection

            % TODO
            mixed = 0;
            vol = 1.0;

        end
        
        function [mixed, vol] = gconic_rec(obj, ifinquiry, ifcyl, dim, r0, r1, ctr0, ctr1, xl, dx)
            % Function: gconic_rec
            % Inputs:
            %   ifinquiry - Inquiry flag to check volume fraction or mixed/clean status
            %   ifcyl - Indicates which axis defines the cylinder (1 = x, 2 = y, 3 = z)
            %   dim - Dimensionality of the problem
            %   r0 - Radius of the cylinder at the base
            %   r1 - Radius of the cylinder at the top
            %   ctr0 - Center of the base of the cylinder
            %   ctr1 - Center of the top of the cylinder
            %   xl - Lower coordinates of the cell
            %   dx - Cell dimensions
            %
            % Outputs:
            %   mixed - Indicates if the cell intersects the cylindrical region (0 = no, 1 = yes)
            %   vol - Volume of the intersection

            % TODO
            mixed = 0;
            vol = 1.0;

        end
        
        function [nnode_intface, nodelist_intface, nnode_out, ...
                  coords1d_out, nface_for_zone_out, nnode_for_face_out, ...
                  nodelist_for_face_out] = poly3d_cut_by_plane(obj, dim, ...
                  norm_direction, plane_location, x0_scale, dx_scale, ...
                  nface, nedge, nnode, coords1d, num_node_for_face, ...
                  nodelist_for_face, edgelist_for_face, nodelist_for_edge)
            % Function: poly3d_cut_by_plane
            % Inputs:
            %   dim - Dimensionality of the problem
            %   norm_direction - Normal direction of the plane
            %   plane_location - Location of the plane
            %   x0_scale - Scaling factor for x0 coordinates
            %   dx_scale - Scaling factor for dx
            %   nface - Number of faces in the polyhedron
            %   nedge - Number of edges in the polyhedron
            %   nnode - Number of nodes in the polyhedron
            %   coords1d - Coordinates of the polyhedron vertices
            %   num_node_for_face - Number of nodes per face
            %   nodelist_for_face - Node indices for each face
            %   edgelist_for_face - Edge indices for each face
            %   nodelist_for_edge - Node indices for each edge
            %
            % Outputs:
            %   nnode_intface - Number of nodes in the interface
            %   nodelist_intface - Node indices for the interface
            %   nnode_out - Number of nodes in the output polyhedron
            %   coords1d_out - Coordinates of the output polyhedron vertices
            %   nface_for_zone_out - Number of faces in the output zone
            %   nnode_for_face_out - Number of nodes per face in the output polyhedron
            %   nodelist_for_face_out - Node indices for each face in the output polyhedron

            % TODO

        end

    end
end
