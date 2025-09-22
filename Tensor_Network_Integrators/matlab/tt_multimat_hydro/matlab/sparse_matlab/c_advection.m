classdef c_advection < handle
    % c_advection   advection methods

    properties
        small = 1e-8;

        % Intercell edge variables
        nmat_for_edge = [];
        matid_for_edge = [];
        dvol_for_edge = [];
        dmass_for_edge = [];
        dener_for_edge = [];

        % Conserved quantities (4D  array, size [nz, ny, nx, nmat])
        vol_for_cell  = [];
        mass_for_cell = [];
        ener_for_cell = [];

    end
    methods

%%%%%%%%%%%%%%% implemented but not tested
        function vel_face = get_face_velocity(~, mesh, dir)
            % get_face_velocity Computes velocity at the cell faces.
            %
            % This method computes the face-centered velocities for 2D or 3D meshes
            % based on the node velocities, depending on the specified direction.
            %
            % Inputs:
            %   mesh - (c_mesh) The mesh object containing grid and velocity data.
            %   dir  - (int) The direction for which the face velocity is computed (0, 1, or 2).
            %
            % Outputs:
            %   vel_face - (array) The computed face-centered velocities.
            %              For 2D, it's a 2D array [ny, nx].
            %              For 3D, it's a 3D array [nz, ny, nx].
            %
            % Uses:
            %   mesh.dim_prob         - (int) Dimensionality of the mesh (2D or 3D).
            %   mesh.ncell_prob       - (array of int) Number of cells in each direction.
            %   mesh.nbdry_prob       - (int) Number of boundary cells.
            %   mesh.vav_for_2dnode   - (3D array of double) Averaged node velocities for 2D meshes.
            %   mesh.vav_for_3dnode   - (4D array of double) Averaged node velocities for 3D meshes.
            %
            % Modifies:
            %   None. This function returns the face velocities without modifying the input data.
            %
            % Original C declaration:
            % void get_face_velocity(int dim, int *ncell, int nbdry, int dir,
            %                        double ***vel_for_2dnode, double ****vel_for_3dnode,
            %                        double ***vel_face_2d, double ****vel_face_3d)

            dim = mesh.dim_prob;  % Dimensionality of the mesh
            ncell = mesh.ncell_prob;
            nbdry = mesh.nbdry_prob;

            % Compute extended cell and node sizes
            ncell_ext = ncell + 2 * nbdry;
            nnode_ext = ncell_ext + 1;

            if dim == 2
                % 2D case
                sizes = [0, 0];
                if dir == 1
                    sizes = [nnode_ext(1), ncell_ext(2)];
                elseif dir == 2
                    sizes = [ncell_ext(1), nnode_ext(2)];
                end
                vel_face = zeros(sizes(2), sizes(1));

                if dir == 1
                    for j = 1:ncell_ext(2)
                        for i = 1:nnode_ext(1)
                            vel_face(j, i) = 0.5 * (mesh.vav_for_2dnode(j, i, dir) + mesh.vav_for_2dnode(j + 1, i, dir));
                        end
                    end
                elseif dir == 2
                    for j = 1:nnode_ext(2)
                        for i = 1:ncell_ext(1)
                            vel_face(j, i) = 0.5 * (mesh.vav_for_2dnode(j, i, dir) + mesh.vav_for_2dnode(j, i + 1, dir));
                        end
                    end
                end

            elseif dim == 3
                % 3D case
                sizes = [0, 0, 0];
                if dir == 0
                    sizes = [nnode_ext(1), ncell_ext(2), ncell_ext(3)];
                elseif dir == 1
                    sizes = [ncell_ext(1), nnode_ext(2), ncell_ext(3)];
                elseif dir == 2
                    sizes = [ncell_ext(1), ncell_ext(2), nnode_ext(3)];
                end
                vel_face = zeros(sizes(3), sizes(2), sizes(1));

                if dir == 0
                    for k = 1:ncell_ext(3)
                        for j = 1:ncell_ext(2)
                            for i = 1:nnode_ext(1)
                                vel_face(k, j, i) = 0.25 * (mesh.vav_for_3dnode(k, j, i, dir) + ...
                                                            mesh.vav_for_3dnode(k, j + 1, i, dir) + ...
                                                            mesh.vav_for_3dnode(k + 1, j, i, dir) + ...
                                                            mesh.vav_for_3dnode(k + 1, j + 1, i, dir));
                            end
                        end
                    end
                elseif dir == 1
                    for k = 1:ncell_ext(3)
                        for j = 1:nnode_ext(2)
                            for i = 1:ncell_ext(1)
                                vel_face(k, j, i) = 0.25 * (mesh.vav_for_3dnode(k, j, i, dir) + ...
                                                            mesh.vav_for_3dnode(k, j, i + 1, dir) + ...
                                                            mesh.vav_for_3dnode(k + 1, j, i, dir) + ...
                                                            mesh.vav_for_3dnode(k + 1, j, i + 1, dir));
                            end
                        end
                    end
                elseif dir == 2
                    for k = 1:nnode_ext(3)
                        for j = 1:ncell_ext(2)
                            for i = 1:ncell_ext(1)
                                vel_face(k, j, i) = 0.25 * (mesh.vav_for_3dnode(k, j, i, dir) + ...
                                                            mesh.vav_for_3dnode(k, j, i + 1, dir) + ...
                                                            mesh.vav_for_3dnode(k, j + 1, i, dir) + ...
                                                            mesh.vav_for_3dnode(k, j + 1, i + 1, dir));
                            end
                        end
                    end
                end
            end
        end

        function courant_adv = advection(obj, mesh, mat, ncycle, dir, dt, ...
                                         btype_lower, btype_upper, is_solid, gamma_ea_mat)
            % advection Performs advection operations on the mesh.
            %
            % This method handles advection operations for 2D and 3D meshes, including
            % face velocity calculations, material data mapping, and boundary conditions.
            %
            % Inputs:
            %   obj          - (c_advection) The current advection object instance.
            %   mesh         - (c_mesh) The mesh object containing grid and material data.
            %   mat          - (c_mat) The material object containing material information.
            %   ncycle       - (int) The current cycle number.
            %   dir          - (int) The direction of advection.
            %   dt           - (double) The time step.
            %   btype_lower  - (array of int) Boundary type for the lower boundaries.
            %   btype_upper  - (array of int) Boundary type for the upper boundaries.
            %   is_solid     - (array of int) Indicator array for solid materials.
            %   gamma_ea_mat - (array of double) Gamma values for each material.
            %
            % Outputs:
            %   courant_adv  - (double) The Courant number for the advection operation.
            %
            % Uses:
            %   mesh.dim_prob        - (int) Dimensionality of the mesh (2D or 3D).
            %   mesh.ncell_prob      - (array of int) Number of cells in each spatial dimension.
            %   mesh.nbdry_prob      - (int) Number of boundary cells.
            %   mesh.xl_prob         - (array of double) Lower bounds of the grid.
            %   mesh.dx_prob         - (array of double) Grid cell sizes.
            %   mat.nmat_prob        - (int) Number of materials in the mesh.
            %   mat.matid_ea_mat     - (array of int) Material IDs for each material.
            %   mat.ijk_in_mixcell   - (cell array) Indices of mixed cells in the mesh.
            %
            % Modifies:
            %   mesh.rho_for_2dcell  - (2D array of double) Updates densities for 2D cells.
            %   mesh.ei_for_2dcell   - (2D array of double) Updates internal energies for 2D cells.
            %   mesh.pres_for_2dcell - (2D array of double) Updates pressures for 2D cells.
            %   mesh.vel_for_2dnode  - (3D array of double) Updates velocities at 2D nodes.
            %   mesh.rho_for_3dcell  - (3D array of double) Updates densities for 3D cells.
            %   mesh.ei_for_3dcell   - (3D array of double) Updates internal energies for 3D cells.
            %   mesh.pres_for_3dcell - (3D array of double) Updates pressures for 3D cells.
            %   mesh.vel_for_3dnode  - (4D array of double) Updates velocities at 3D nodes.
            %   mesh.cs_for_2dcell   - (2D array of double) Updates sound speeds for 2D cells.
            %   mesh.cs_for_3dcell   - (3D array of double) Updates sound speeds for 3D cells.
            %
            % Original C declaration:
            % void advection(int fileid, int dim, int *ncell, int nbdry, double *xl_prob, double *dx_prob,
            %                int nmat_prob, int *matid_ea_mat, int *is_solid, double *gamma_ea_mat,
            %                int ncycle, int dir, double dt,
            %                Bdry_Type *btype_lower, Bdry_Type *btype_upper,
            %                double **rho_for_2dcell, double **ei_for_2dcell, double **pres_for_2dcell,
            %                double ***vel_for_2dnode, double ***vav_for_2dnode,
            %                int  **nmat_for_2dcell, int  **matid_for_2dcell,
            %                double ***rho_for_3dcell, double ***ei_for_3dcell, double ***pres_for_3dcell,
            %                double ****vel_for_3dnode, double ****vav_for_3dnode,
            %                int ***nmat_for_3dcell, int ***matid_for_3dcell,
            %                double *courant_adv)
            global debug xio;

            % Get mesh properties
            dim = mesh.dim_prob;
            ncell = mesh.ncell_prob;
            nbdry = mesh.nbdry_prob;

            % Initialize variables
            ncell_ext = ncell + 2 * nbdry;  % Extended cell size (including boundaries)

            % Get face velocity (single output)
            vel_face = obj.get_face_velocity(mesh, dir);
            if dim == 2
                vel_face_2d = vel_face;
                vel_face_3d = [];
            else
                vel_face_2d = [];
                vel_face_3d = vel_face;
            end

            % 2D case
            if dim == 2
                mesh.mixcell_for_2dcell = -ones(ncell_ext(2), ncell_ext(1));
                for mx = 1:mat.nmixcell
                    i = mat.ijk_in_mixcell(mx,1);
                    j = mat.ijk_in_mixcell(mx,2);
                    mesh.mixcell_for_2dcell(j, i) = mx;
                end
                mesh.mixcell_for_3dcell = [];

            % 3D case
            elseif dim == 3
                mesh.mixcell_for_3dcell = -ones(ncell_ext(3), ncell_ext(2), ncell_ext(1));
                for mx = 1:mat.nmixcell
                    i = mat.ijk_in_mixcell(mx,1);
                    j = mat.ijk_in_mixcell(mx,2);
                    k = mat.ijk_in_mixcell(mx,3);
                    mesh.mixcell_for_3dcell(k, j, i) = mx;
                end
                mesh.mixcell_for_2dcell = [];
            end

            % Get material polygons/polyhedrons
            mat.get_mpoly(mesh);

%if (debug), xio.write_dump(mesh, mat, 37707, 1.0), end
            % Perform material mapping and update the Courant number
            courant_adv = obj.mapping(mesh, mat, dir, dt, ...
                          btype_lower, btype_upper, is_solid, gamma_ea_mat, ...
                          vel_face_2d, vel_face_3d);

%if (debug)
%    xio.write_dump(mesh, mat, 37708, 1.1);
%end
        end

        function courant_adv = mapping(obj, mesh, mat, dir, dt, ...
                                       btype_lower, btype_upper, is_solid, gamma_ea_mat, ...
                                       vel_face_2d, vel_face_3d)
            % mapping Switches between 2D and 3D material mapping for advection.
            %
            % This method switches between 2D and 3D mapping functions based on the
            % dimensionality of the mesh. It handles the advection and material mapping
            % operations, updating the state variables and calculating the Courant number.
            %
            % Inputs:
            %   obj           - (c_advection) The current advection object instance.
            %   mesh          - (c_mesh) The mesh object containing grid and material data.
            %   mat           - (c_mat) The material object containing material information.
            %   dir           - (int) The direction of advection.
            %   dt            - (double) The time step.
            %   btype_lower   - (array of int) Boundary type for the lower boundaries.
            %   btype_upper   - (array of int) Boundary type for the upper boundaries.
            %   is_solid      - (array of int) Indicator array for solid materials.
            %   gamma_ea_mat  - (array of double) Gamma values for each material.
            %   vel_face_2d   - (2D array of double) Face velocities for 2D cells.
            %   vel_face_3d   - (3D array of double) Face velocities for 3D cells.
            %
            % Outputs:
            %   courant_adv   - (double) The Courant number for the advection operation.
            %
            % Uses:
            %   mesh.dim_prob         - (int) Dimensionality of the mesh (2D or 3D).
            %   mesh.ncell_prob       - (array of int) Number of cells in each spatial dimension.
            %   mesh.nbdry_prob       - (int) Number of boundary cells.
            %   mesh.xl_prob          - (array of double) Lower bounds of the grid.
            %   mesh.dx_prob          - (array of double) Grid cell sizes.
            %   mat.nmat_prob         - (int) Number of materials in the mesh.
            %   mat.matid_ea_mat      - (array of int) Material IDs for each material.
            %
            % Modifies:
            %   mesh.rho_for_2dcell   - (2D array of double) Updates densities for 2D cells.
            %   mesh.ei_for_2dcell    - (2D array of double) Updates internal energies for 2D cells.
            %   mesh.pres_for_2dcell  - (2D array of double) Updates pressures for 2D cells.
            %   mesh.vel_for_2dnode   - (3D array of double) Updates velocities at 2D nodes.
            %   mesh.rho_for_3dcell   - (3D array of double) Updates densities for 3D cells.
            %   mesh.ei_for_3dcell    - (3D array of double) Updates internal energies for 3D cells.
            %   mesh.pres_for_3dcell  - (3D array of double) Updates pressures for 3D cells.
            %   mesh.vel_for_3dnode   - (4D array of double) Updates velocities at 3D nodes.
            %
            % Original C declaration:
            % void mapping(int dim, int *ncell, int nbdry, double *xl_prob, double *dx_prob, int dir, double dt,
            %              Bdry_Type *btype_lower, Bdry_Type *btype_upper,
            %              int nmat_prob, int *matid_ea_mat, int *is_solid, double *gamma_ea_mat,
            %              double **rho_for_2dcell, double **ei_for_2dcell, double **pres_for_2dcell,
            %              double ***vel_for_2dnode, double **vel_face_2d,
            %              int  **nmat_for_2dcell, int  **matid_for_2dcell, int  **mixcell_for_2dcell,
            %              double ***rho_for_3dcell, double ***ei_for_3dcell, double ***pres_for_3dcell,
            %              double ****vel_for_3dnode, double ***vel_face_3d,
            %              int ***nmat_for_3dcell, int ***matid_for_3dcell, int ***mixcell_for_3dcell,
            %              double *courant_adv)

            % Get mesh properties
            dim = mesh.dim_prob;

            % Call mapping2d or mapping3d based on dimensionality
            if dim == 2
                % Call 2D mapping
                courant_adv = obj.mapping2d(mesh, mat, dir, dt, ...
                                            btype_lower, btype_upper, is_solid, ...
                                            gamma_ea_mat, vel_face_2d);
            elseif dim == 3
                % Call 3D mapping
                courant_adv = obj.mapping3d(mesh, mat, dir, dt, ...
                                            btype_lower, btype_upper, is_solid, ...
                                            gamma_ea_mat, vel_face_3d);
            end
        end

        function courant_adv = mapping2d(~, mesh, mat, dir, dt, ...
                                         btype_lower, btype_upper, is_solid, gamma_ea_mat, ...
                                         vel_face_2d)
            % mapping2d Performs material mapping and advection for 2D meshes.
            %
            % This method handles the material mapping and advection operations for
            % 2D meshes. It updates the cell velocities, densities, pressures, and
            % internal energies, while calculating the Courant number for stability.
            %
            % Inputs:
            %   obj             - (c_advection) The current advection object instance.
            %   mesh            - (c_mesh) The mesh object containing grid and material data.
            %   mat             - (c_mat) The material object containing material information.
            %   dir             - (int) The direction of advection.
            %   dt              - (double) The time step.
            %   btype_lower     - (array of int) Boundary type for the lower boundaries.
            %   btype_upper     - (array of int) Boundary type for the upper boundaries.
            %   is_solid        - (array of int) Indicator array for solid materials.
            %   gamma_ea_mat    - (array of double) Gamma values for each material.
            %   vel_face_2d     - (2D array of double) Face velocities for 2D cells.
            %
            % Outputs:
            %   courant_adv     - (double) The Courant number for the advection operation.
            %
            % Uses:
            %   mesh.ncell_prob      - (array of int) Number of cells in each spatial dimension.
            %   mesh.nbdry_prob      - (int) Number of boundary cells.
            %   mesh.xl_prob         - (array of double) Lower bounds of the grid.
            %   mesh.dx_prob         - (array of double) Grid cell sizes.
            %   mat.nmat_prob        - (int) Number of materials in the mesh.
            %   mat.matid_ea_mat     - (array of int) Material IDs for each material.
            %
            % Modifies:
            %   mesh.rho_for_2dcell  - (2D array of double) Updates densities for 2D cells.
            %   mesh.ei_for_2dcell   - (2D array of double) Updates internal energies for 2D cells.
            %   mesh.pres_for_2dcell - (2D array of double) Updates pressures for 2D cells.
            %   mesh.vel_for_2dnode  - (3D array of double) Updates velocities at 2D nodes.
            %
            % Original C declaration:
            % void mapping2d(int *ncell, int nbdry, double *xl_prob, double *dx_prob, int dir, double dt,
            %                Bdry_Type *btype_lower, Bdry_Type *btype_upper, int nmat_prob, int *matid_ea_mat,
            %                int *is_solid, double *gamma_ea_mat, double **vel_face_2d,
            %                double **rho_for_2dcell, double **ei_for_2dcell, double **pres_for_2dcell,
            %                double ***vel_for_2dnode, int **nmat_for_2dcell, int **matid_for_2dcell,
            %                int **mixcell_for_2dcell, double *courant_adv)

            global bdry small eos;
            small = 1e-8;

            % Initialize arrays and scalar variables
            j_start = zeros(1, 2);  % Equivalent to int j_start[2] in C
            j_end = zeros(1, 2);    % Equivalent to int j_end[2] in C
            i_start = zeros(1, 2);  % Equivalent to int i_start[2] in C
            i_end = zeros(1, 2);    % Equivalent to int i_end[2] in C

            % Logical flags (can also be doubles if needed)
            if_mom_conserved = true;    % Logical variable equivalent to int with value 1
            if_mom_homoginized = true;  % Logical variable equivalent to int with value 1


            % Get mesh properties
            ncell = mesh.ncell_prob;
            nbdry = mesh.nbdry_prob;
            xl = mesh.xl_prob;
            dx = mesh.dx;
            vol_cell = prod(dx);
            dim = 2;

            ncell_ext = ncell + 2 * nbdry;
            ncell_bdry = ncell_ext - nbdry;    % cell at the right boundary
            nnode_ext = ncell_ext + 1;
            nnode_bdry = nnode_ext - nbdry;    % node at the right boundary
            hdx = 0.5 * dx;
            sizes_edge = ncell_ext;                % Size of edges in the x and y directions
            sizes_edge(dir) = nnode_ext(dir);      % Adjust size for the direction of advection

            % Material properties
            nmat_prob = mat.nmat_prob;
            matid_ea_mat = mat.matids_prob;

            % Allocate memory for necessary variables
            matid_adv = zeros(2*nmat_prob, 1);
            indice_adv = zeros(nmat_prob, 1);  % Pointer to the second half of matid_adv
            vol_adv = zeros(nmat_prob, 1);
            mass_adv = zeros(nmat_prob, 1);
            ener_adv = zeros(nmat_prob, 1);

            % Compute average density
            rho_average = mean(mesh.rho_for_2dcell(nbdry+1:ncell_bdry(2), nbdry+1:ncell_bdry(1)), 'all');

            %   calculate the vol, mass, and energy advected across the cell interfaces.
            nmat_for_edge = zeros(sizes_edge(2), sizes_edge(1));  % Initialize the material count at each edge to 0

            lsize_max = sizes_edge(1) * sizes_edge(2) * 2;  % Maximum size of materials crossing the edges
            dvol_for_edge  = zeros(lsize_max, 1);  % Array for volumes crossing edges
            dmass_for_edge = zeros(lsize_max, 1);  % Array for masses crossing edges
            dener_for_edge = zeros(lsize_max, 1);  % Array for energies crossing edges
            matid_for_edge = zeros(lsize_max, 1);  % Array for material IDs crossing edges

            loc_edge = 1;  % Initialize current location in dmass_for_edge, etc.
            courant_adv = 0.0; % Initialize courant number and velocity arrays

            if (dir == 1)
                for j = nbdry+1:ncell_bdry(2)
                    jc = j;
                    xl_cell(2) = xl(2) + (j - nbdry - 1) * dx(2);
                    xl_slab(2) = xl_cell(2);  % for the slab advected
                    xr_slab(2) = xl_cell(2) + dx(2);
                    inward_norm(2) = 0.0;

                    for edge = nbdry+1:nnode_bdry(1)
                        x_edge = xl(dir) + (edge - nbdry - 1) * dx(dir);

                        dist = vel_face_2d(j, edge) * dt;
                        frac = dist / dx(dir);
                        courant_adv = max(courant_adv, abs(frac));
                        if courant_adv > 1
                            warning("courant_adv = %e in mapping2d", courant_adv);
                        end
                        clower = edge - 1;
                        cupper = edge;

                        if frac > small  % cell clower moves the following mass and energy to cupper
                            if mesh.nmat_for_2dcell(j, clower) == 1
                                % Reallocate arrays if necessary
                                if loc_edge + 1 > lsize_max
                                    lsize_max = round(1.2 * lsize_max) + 8;
                                    dvol_for_edge = [dvol_for_edge; zeros(lsize_max - numel(dvol_for_edge), 1)];
                                    dmass_for_edge = [dmass_for_edge; zeros(lsize_max - numel(dmass_for_edge), 1)];
                                    dener_for_edge = [dener_for_edge; zeros(lsize_max - numel(dener_for_edge), 1)];
                                end
                                nmat_for_edge(j, edge) = 1;

                                dvol = dx(2) * dist;
                                dmass = dvol * mesh.rho_for_2dcell(j, clower);
                                dener = dvol * mesh.ei_for_2dcell(j, clower);

                                matid_for_edge(loc_edge) = mesh.matid_for_2dcell(j, clower);
                                dvol_for_edge(loc_edge) = dvol;
                                dmass_for_edge(loc_edge) = dmass;
                                dener_for_edge(loc_edge) = dener;
                                loc_edge = loc_edge + 1;
                            else
                                % Mixed materials case
                                ic = clower;
                                xl_slab(dir) = x_edge - abs(dist);
                                xr_slab(dir) = x_edge;
                                xl_cell(dir) = x_edge - dx(1);
                                inward_norm(dir) = 1.0;  % inward norm within the slab

                                slab_faceid = 0;  % xl-face of slab
                                mixcell = mesh.mixcell_for_2dcell(jc, ic);
                                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] ...
                                    = mat.advect2d(xl_cell, dx, mixcell, xl_slab, ...
                                                   xr_slab, inward_norm, slab_faceid);

                                if loc_edge + nmat_adv > lsize_max
                                    lsize_max = round(1.2 * lsize_max) + nmat_adv;
                                    dvol_for_edge = [dvol_for_edge; zeros(lsize_max - numel(dvol_for_edge), 1)];
                                    dmass_for_edge = [dmass_for_edge; zeros(lsize_max - numel(dmass_for_edge), 1)];
                                    dener_for_edge = [dener_for_edge; zeros(lsize_max - numel(dener_for_edge), 1)];
                                end
                                nmat_for_edge(j, edge) = nmat_adv;
                                matid_for_edge(loc_edge:loc_edge+nmat_adv-1) = matid_adv(1:nmat_adv);
                                dvol_for_edge(loc_edge:loc_edge+nmat_adv-1) = vol_adv(1:nmat_adv);
                                dmass_for_edge(loc_edge:loc_edge+nmat_adv-1) = mass_adv(1:nmat_adv);
                                dener_for_edge(loc_edge:loc_edge+nmat_adv-1) = ener_adv(1:nmat_adv);
                                loc_edge = loc_edge + nmat_adv;
                            end
                        elseif frac < -small
                            adist = abs(dist);
                            if mesh.nmat_for_2dcell(j, cupper) == 1
                                % Reallocate arrays if necessary
                                if loc_edge + 1 > lsize_max
                                    lsize_max = round(1.2 * lsize_max) + 8;
                                    dvol_for_edge = [dvol_for_edge; zeros(lsize_max - numel(dvol_for_edge), 1)];
                                    dmass_for_edge = [dmass_for_edge; zeros(lsize_max - numel(dmass_for_edge), 1)];
                                    dener_for_edge = [dener_for_edge; zeros(lsize_max - numel(dener_for_edge), 1)];
                                end
                                nmat_for_edge(j, edge) = 1;
                                dvol = dx(2) * adist;
                                dmass = dvol * mesh.rho_for_2dcell(j, cupper);
                                dener = dvol * mesh.ei_for_2dcell(j, cupper);

                                matid_for_edge(loc_edge) = mesh.matid_for_2dcell(j, cupper);
                                dvol_for_edge(loc_edge) = -dvol;
                                dmass_for_edge(loc_edge) = -dmass;
                                dener_for_edge(loc_edge) = -dener;
                                loc_edge = loc_edge + 1;
                            else
                                % Mixed materials case for negative frac
                                ic = cupper;
                                xl_slab(dir) = x_edge;
                                xr_slab(dir) = x_edge + abs(dist);
                                xl_cell(dir) = x_edge;
                                inward_norm(dir) = -1.0;  % inward norm within the slab

                                slab_faceid = 1;  % xr-face of slab
                                mixcell = mesh.mixcell_for_2dcell(jc, ic);
                                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] ...
                                    = mat.advect2d(xl_cell, dx, mixcell, xl_slab, ...
                                                   xr_slab, inward_norm, slab_faceid);

                                if loc_edge + nmat_adv > lsize_max
                                    lsize_max = round(1.2 * lsize_max) + nmat_adv;
                                    dvol_for_edge = [dvol_for_edge; zeros(lsize_max - numel(dvol_for_edge), 1)];
                                    dmass_for_edge = [dmass_for_edge; zeros(lsize_max - numel(dmass_for_edge), 1)];
                                    dener_for_edge = [dener_for_edge; zeros(lsize_max - numel(dener_for_edge), 1)];
                                end
                                nmat_for_edge(j, edge) = nmat_adv;
                                for idx = 1:nmat_adv
                                    matid_for_edge(loc_edge) = matid_adv(idx);
                                    dvol_for_edge(loc_edge) = -vol_adv(idx);
                                    dmass_for_edge(loc_edge) = -mass_adv(idx);
                                    dener_for_edge(loc_edge) = -ener_adv(idx);
                                    loc_edge = loc_edge + 1;
                                end
                            end
                        end
                    end
                end

            elseif (dir == 2)
                % Loop over edges and calculate advection
                for edge = nbdry+1:nnode_bdry(2)
                    y_edge = xl(dir) + (edge - nbdry) * dx(dir);  % y-coordinate of the edge

                    for i = nbdry+1:ncell_bdry(1)
                        ic = i;
                        xl_cell(1) = xl(1) + (i - nbdry) * dx(1);
                        xl_slab(1) = xl_cell(1);  % for the slab being advected
                        xr_slab(1) = xl_cell(1) + dx(1);

                        inward_norm = zeros(1, 2);  % initialize inward norm

                        dist = vel_face_2d(edge, i) * dt;
                        frac = dist / dx(dir);
                        courant_adv = max(courant_adv, abs(frac));
                        assert(courant_adv < 1.0);  % Ensure Courant number is less than 1

                        clower = edge - 1;
                        cupper = edge;

                        if frac > small  % cell clower moves mass and energy to cupper
                            if mesh.nmat_for_2dcell(clower, i) == 1
                                % Increase array size dynamically if needed
                                if loc_edge + 1 > lsize_max
                                    lsize_max = round(1.2 * lsize_max) + 8;
                                    dvol_for_edge(lsize_max) = 0;  % Expand array in MATLAB
                                    dmass_for_edge(lsize_max) = 0;
                                    dener_for_edge(lsize_max) = 0;
                                end
                                nmat_for_edge(edge, i) = 1;
                                dvol = dx(1) * dist;
                                dmass = dvol * mesh.rho_for_2dcell(clower, i);
                                dener = dvol * mesh.ei_for_2dcell(clower, i);

                                matid_for_edge(loc_edge) = mesh.matid_for_2dcell(clower, i);
                                dvol_for_edge(loc_edge) = dvol;
                                dmass_for_edge(loc_edge) = dmass;
                                dener_for_edge(loc_edge) = dener;
                                loc_edge = loc_edge + 1;
                            else
                                % Handle mixed cell case
                                jc = clower;
                                xl_slab(dir) = y_edge - abs(dist);
                                xr_slab(dir) = y_edge;
                                xl_cell(dir) = y_edge - dx(dir);
                                inward_norm(dir) = 1.0;  % inward norm within the slab

                                slab_faceid = 2;  % yl-face of slab
                                mixcell = mesh.mixcell_for_2dcell(jc, ic);
                                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] ...
                                    = mat.advect2d(xl_cell, dx, mixcell, xl_slab, ...
                                                   xr_slab, inward_norm, slab_faceid);

                                if loc_edge + nmat_adv > lsize_max
                                    lsize_max = round(1.2 * lsize_max) + nmat_adv;
                                    dvol_for_edge(lsize_max) = 0;  % Expand array in MATLAB
                                    dmass_for_edge(lsize_max) = 0;
                                    dener_for_edge(lsize_max) = 0;
                                end
                                nmat_for_edge(edge, i) = nmat_adv;
                                matid_for_edge(loc_edge:loc_edge + nmat_adv - 1) = matid_adv;
                                dvol_for_edge(loc_edge:loc_edge + nmat_adv - 1) = vol_adv;
                                dmass_for_edge(loc_edge:loc_edge + nmat_adv - 1) = mass_adv;
                                dener_for_edge(loc_edge:loc_edge + nmat_adv - 1) = ener_adv;
                                loc_edge = loc_edge + nmat_adv;
                            end
                        elseif frac < -small  % cell cupper moves mass and energy to clower
                            adist = abs(dist);
                            if mesh.nmat_for_2dcell(cupper, i) == 1
                                if loc_edge + 1 > lsize_max
                                    lsize_max = round(1.2 * lsize_max) + 8;
                                    dvol_for_edge(lsize_max) = 0;  % Expand array in MATLAB
                                    dmass_for_edge(lsize_max) = 0;
                                    dener_for_edge(lsize_max) = 0;
                                end
                                nmat_for_edge(edge, i) = 1;
                                dvol = dx(1) * adist;
                                dmass = dvol * mesh.rho_for_2dcell(cupper, i);
                                dener = dvol * mesh.ei_for_2dcell(cupper, i);

                                matid_for_edge(loc_edge) = mesh.matid_for_2dcell(cupper, i);
                                dvol_for_edge(loc_edge) = -dvol;
                                dmass_for_edge(loc_edge) = -dmass;
                                dener_for_edge(loc_edge) = -dener;
                                loc_edge = loc_edge + 1;
                            else
                                % Handle mixed cell case for cupper
                                jc = cupper;
                                xl_slab(dir) = y_edge;
                                xr_slab(dir) = y_edge + abs(dist);
                                xl_cell(dir) = y_edge;
                                inward_norm(dir) = -1.0;  % inward norm within the slab

                                slab_faceid = 3;  % yr-face of slab
                                mixcell = mesh.mixcell_for_2dcell(jc, ic);
                                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] ...
                                    = mat.advect2d(xl_cell, dx, mixcell, xl_slab, ...
                                                   xr_slab, inward_norm, slab_faceid);

                                if loc_edge + nmat_adv > lsize_max
                                    lsize_max = round(1.2 * lsize_max) + nmat_adv;
                                    dvol_for_edge(lsize_max) = 0;  % Expand array in MATLAB
                                    dmass_for_edge(lsize_max) = 0;
                                    dener_for_edge(lsize_max) = 0;
                                end
                                nmat_for_edge(edge, i) = nmat_adv;
                                for idx = 1:nmat_adv
                                    matid_for_edge(loc_edge) = matid_adv(idx);
                                    dvol_for_edge(loc_edge) = -vol_adv(idx);
                                    dmass_for_edge(loc_edge) = -mass_adv(idx);
                                    dener_for_edge(loc_edge) = -ener_adv(idx);
                                    loc_edge = loc_edge + 1;
                                end
                            end
                        end
                    end  % end i-loop
                end  % end edge-loop

            end

            % Allocate 2D array for loc_for_2dedge
            loc_for_2dedge = -ones(sizes_edge(2), sizes_edge(1));  % Initialize with -1 to handle zero-material cases
            loc_edge = 1;  % Initialize the current offset (may need to be set to 1)

            % Loop over the edges and update loc_for_2dedge
            for j = 1:sizes_edge(2)
                for i = 1:sizes_edge(1)
                    if nmat_for_edge(j, i) > 0
                        loc_for_2dedge(j, i) = loc_edge;
                        loc_edge = loc_edge + nmat_for_edge(j, i);
                    end
                end
            end

            % Initialize material count for cells
            bdry.nmat_for_cell = zeros(ncell_ext(2), ncell_ext(1));
            bdry.nmat_for_cell(:,:) = mesh.nmat_for_2dcell(:,:);

            % Find the max number of materials in each cell (x-direction)
            if dir == 1
                for j = 1:sizes_edge(2)
                    jc = j;
                    for i = 1:sizes_edge(1)
                        nm_cross_edge = nmat_for_edge(j, i);
                        loc_edge = loc_for_2dedge(j, i);
                        if loc_edge < 0, continue; end

                        if dvol_for_edge(loc_edge) > 0.0
                            ic = i;
                        else
                            ic = i - 1;
                        end

                        nm_cell = mesh.nmat_for_2dcell(jc, ic);
                        if nm_cell == 1
                            matids = mesh.matid_for_2dcell(jc, ic);
                        else
                            mx = mesh.mixcell_for_2dcell(jc, ic);
                            matids = mat.matids_in_mixcell{mx};
                        end

                        for idx = 1:nm_cross_edge
                            m = matid_for_edge(loc_edge);
                            found = false;
                            for s = 1:nm_cell
                                if matids(s) == m
                                    found = true;
                                    break;
                                end
                            end
                            if ~found
                                bdry.nmat_for_cell(jc, ic) = bdry.nmat_for_cell(jc, ic) + 1;
                            end
                            loc_edge = loc_edge + 1;
                        end
                    end
                end


            elseif dir == 2
                % Find the max number of materials in each cell (y-direction)
                for j = 1:sizes_edge(2)
                    for i = 1:sizes_edge(1)
                        ic = i;

                        nm_cross_edge = nmat_for_edge(j, i);
                        loc_edge = loc_for_2dedge(j, i);
                        if loc_edge < 0, continue; end

                        if dvol_for_edge(loc_edge) > 0.0
                            jc = j;
                        else
                            jc = j - 1;
                        end

                        nm_cell = mesh.nmat_for_2dcell(jc, ic);
                        if nm_cell == 1
                            matids = mesh.matid_for_2dcell(jc, ic);
                        else
                            mx = mesh.mixcell_for_2dcell(jc, ic);
                            matids = mat.matids_in_mixcell{mx};
                        end

                        for idx = 1:nm_cross_edge
                            m = matid_for_edge(loc_edge);
                            found = false;
                            for s = 1:nm_cell
                                if matids(s) == m
                                    found = true;
                                    break;
                                end
                            end
                            if ~found
                                bdry.nmat_for_cell(jc, ic) = bdry.nmat_for_cell(jc, ic) + 1;
                            end
                            loc_edge = loc_edge + 1;
                        end
                    end
                end
            end

            % apply boundary condition to bdry.nmat_for_cell
            bdry.bdry_cell_1var_2d(mesh, btype_lower, btype_upper);

            bdry.matid_for_cell = cell(ncell_ext(2), ncell_ext(1));  % 2D cell array for material IDs
            bdry.vol_for_cell = cell(ncell_ext(2), ncell_ext(1));    % 2D cell array for volumes
            bdry.mass_for_cell = cell(ncell_ext(2), ncell_ext(1));   % 2D cell array for masses
            bdry.ener_for_cell = cell(ncell_ext(2), ncell_ext(1));   % 2D cell array for energies

            for j = 1:ncell_ext(2)
                for i = 1:ncell_ext(1)
                    nm_cell = bdry.nmat_for_cell(j, i);
                    bdry.matid_for_cell{j, i} = zeros(1, nm_cell);
                    bdry.vol_for_cell{j, i} = zeros(1, nm_cell);
                    bdry.mass_for_cell{j, i} = zeros(1, nm_cell);
                    bdry.ener_for_cell{j, i} = zeros(1, nm_cell);

                    nm_cell = mesh.nmat_for_2dcell(j, i);
                    bdry.nmat_for_cell(j, i) = nm_cell;

                    % If there is only one material in the cell
                    if nm_cell == 1
                        bdry.matid_for_cell{j, i}(1) = mesh.matid_for_2dcell(j, i);
                        bdry.vol_for_cell{j, i}(1) = vol_cell;
                        bdry.mass_for_cell{j, i}(1) = mesh.rho_for_2dcell(j, i) * vol_cell;
                        bdry.ener_for_cell{j, i}(1) = mesh.ei_for_2dcell(j, i) * vol_cell;
                    else
                        % If the cell contains mixed materials
                        mx = mesh.mixcell_for_2dcell(j, i);
                        for idx = 1:nm_cell
                            bdry.matid_for_cell{j, i}(idx) = mat.matids_in_mixcell{mx}(idx);
                            dvol = vol_cell * mat.vf_in_mixcell{mx}(idx);
                            bdry.vol_for_cell{j, i}(idx) = dvol;
                            bdry.mass_for_cell{j, i}(idx) = dvol * mat.rho_in_mixcell{mx}(idx);
                            bdry.ener_for_cell{j, i}(idx) = dvol * mat.ei_in_mixcell{mx}(idx);
                        end
                    end
                end
            end

            %  calculate materials after the advection
            if dir == 1
                for j = 1:sizes_edge(2)
                    jc = j;  % Current edge index in y-direction
                    for i = 1:sizes_edge(1)
                        nm_cross_edge = nmat_for_edge(j, i);
                        loc_edge = loc_for_2dedge(j, i);
                        if loc_edge < 0, continue; end

                        % Determine the cells involved in the material exchange
                        if dvol_for_edge(loc_edge) > 0.0
                            ic_out = i - 1;  % The cell where materials are deducted
                            ic_in  = i;      % The cell where materials are added
                        else
                            ic_out = i;      % The cell where materials are deducted
                            ic_in  = i - 1;  % The cell where materials are added
                        end

                        % Get the number of materials and material IDs in the incoming cell
                        nm_cell = bdry.nmat_for_cell(jc, ic_in);
                        matids = bdry.matid_for_cell{jc, ic_in};

                        nm_new = nm_cell;  % Start with the current number of materials

                        % Loop through the materials crossing the edge
                        for idx = 1:nm_cross_edge
                            m = matid_for_edge(loc_edge);
                            found = false;

                            % Check if the material is already present in the incoming cell
                            for s = 1:nm_cell
                                if matids(s) == m
                                    found = true;
                                    break;
                                end
                            end

                            if ~found
                                % Add new material to the incoming cell
                                nm_new = nm_new + 1;
                                bdry.matid_for_cell{jc, ic_in}(nm_new) = m;
                                bdry.vol_for_cell{jc, ic_in}(nm_new) = abs(dvol_for_edge(loc_edge));
                                bdry.mass_for_cell{jc, ic_in}(nm_new) = abs(dmass_for_edge(loc_edge));
                                bdry.ener_for_cell{jc, ic_in}(nm_new) = abs(dener_for_edge(loc_edge));
                                bdry.nmat_for_cell(jc, ic_in) = bdry.nmat_for_cell(jc, ic_in) + 1;
                            else
                                % Update the material properties in the incoming cell
                                bdry.vol_for_cell{jc, ic_in}(s) = bdry.vol_for_cell{jc, ic_in}(s) + abs(dvol_for_edge(loc_edge));
                                bdry.mass_for_cell{jc, ic_in}(s) = bdry.mass_for_cell{jc, ic_in}(s) + abs(dmass_for_edge(loc_edge));
                                bdry.ener_for_cell{jc, ic_in}(s) = bdry.ener_for_cell{jc, ic_in}(s) + abs(dener_for_edge(loc_edge));
                            end

                            loc_edge = loc_edge + 1;
                        end

                        % Update the outgoing cell (deduct materials)
                        loc_edge = loc_for_2dedge(j, i);
                        nm_cell = bdry.nmat_for_cell(jc, ic_out);
                        matids = bdry.matid_for_cell{jc, ic_out};

                        % Loop through the materials crossing the edge and deduct from outgoing cell
                        for idx = 1:nm_cross_edge
                            m = matid_for_edge(loc_edge);
                            found = false;

                            % Check if the material is present in the outgoing cell
                            for s = 1:nm_cell
                                if matids(s) == m
                                    found = true;
                                    break;
                                end
                            end

                            assert(found, 'Material not found in outgoing cell.');

                            % Deduct the material properties from the outgoing cell
                            bdry.vol_for_cell{jc, ic_out}(s) = bdry.vol_for_cell{jc, ic_out}(s) - abs(dvol_for_edge(loc_edge));
                            bdry.mass_for_cell{jc, ic_out}(s) = bdry.mass_for_cell{jc, ic_out}(s) - abs(dmass_for_edge(loc_edge));
                            bdry.ener_for_cell{jc, ic_out}(s) = bdry.ener_for_cell{jc, ic_out}(s) - abs(dener_for_edge(loc_edge));

                            loc_edge = loc_edge + 1;
                        end
                    end
                end

            elseif dir == 2
                for j = 1:sizes_edge(2)
                    for i = 1:sizes_edge(1)
                        ic = i;

                        % Get the number of materials crossing the edge and the current location
                        nm_cross_edge = nmat_for_edge(j, i);
                        loc_edge = loc_for_2dedge(j, i);
                        if loc_edge < 0, continue; end

                        % Determine which cells to add or deduct materials from
                        if dvol_for_edge(loc_edge) > 0.0
                            jc_out = j - 1;  % the cell where materials are deducted
                            jc_in  = j;      % the cell where materials are added
                        else
                            jc_out = j;      % the cell where materials are deducted
                            jc_in  = j - 1;  % the cell where materials are added
                        end

                        % Get materials for the incoming cell
                        nm_cell = bdry.nmat_for_cell(jc_in, ic);
                        matids = bdry.matid_for_cell{jc_in, ic};

                        nm_new = nm_cell;
                        for idx = 1:nm_cross_edge
                            m = matid_for_edge(loc_edge);
                            found = false;
                            for s = 1:nm_cell
                                if matids(s) == m
                                    found = true;
                                    break;
                                end
                            end
                            if ~found
                                % Add new material to the cell
                                nm_new = nm_new + 1;
                                bdry.matid_for_cell{jc_in, ic}(nm_new) = m;
                                bdry.vol_for_cell{jc_in, ic}(nm_new) = abs(dvol_for_edge(loc_edge));
                                bdry.mass_for_cell{jc_in, ic}(nm_new) = abs(dmass_for_edge(loc_edge));
                                bdry.ener_for_cell{jc_in, ic}(nm_new) = abs(dener_for_edge(loc_edge));
                                bdry.nmat_for_cell(jc_in, ic) = bdry.nmat_for_cell(jc_in, ic) + 1;
                            else
                                % Update existing material in the cell
                                bdry.vol_for_cell{jc_in, ic}(s) = bdry.vol_for_cell{jc_in, ic}(s) + abs(dvol_for_edge(loc_edge));
                                bdry.mass_for_cell{jc_in, ic}(s) = bdry.mass_for_cell{jc_in, ic}(s) + abs(dmass_for_edge(loc_edge));
                                bdry.ener_for_cell{jc_in, ic}(s) = bdry.ener_for_cell{jc_in, ic}(s) + abs(dener_for_edge(loc_edge));
                            end
                            loc_edge = loc_edge + 1;
                        end

                        % Deduct materials from the outgoing cell
                        loc_edge = loc_for_2dedge(j, i);
                        nm_cell = bdry.nmat_for_cell(jc_out, ic);
                        matids = bdry.matid_for_cell{jc_out, ic};

                        for idx = 1:nm_cross_edge
                            m = matid_for_edge(loc_edge);
                            found = false;
                            for s = 1:nm_cell
                                if matids(s) == m
                                    found = true;
                                    break;
                                end
                            end
                            assert(found);

                            % Deduct volume, mass, and energy from the outgoing cell
                            bdry.vol_for_cell{jc_out, ic}(s) = bdry.vol_for_cell{jc_out, ic}(s) - abs(dvol_for_edge(loc_edge));
                            bdry.mass_for_cell{jc_out, ic}(s) = bdry.mass_for_cell{jc_out, ic}(s) - abs(dmass_for_edge(loc_edge));
                            bdry.ener_for_cell{jc_out, ic}(s) = bdry.ener_for_cell{jc_out, ic}(s) - abs(dener_for_edge(loc_edge));

                            loc_edge = loc_edge + 1;
                        end
                    end
                end
            end

            %  take out small volume fraction
            nmixcell_int_new = 0; % new mix cell count
            lsize_mx = 0;     % size of the material list in the mixed cells

            for j = nbdry+1:ncell_bdry(2)
                for i = nbdry+1:ncell_bdry(1)
                    nm_new = bdry.nmat_for_cell(j, i);
                    nm_cell = 0;
                    vol_sum = 0.0;

                    % Loop through materials in the current cell
                    for idx = 1:nm_new
                        frac = bdry.vol_for_cell{j, i}(idx) / vol_cell; % OK: for {j, i} == {3, 53} and idx == 2, frac is zero!

                        % If volume fraction is above the threshold
                        if frac > small
                            nm_cell = nm_cell + 1;
                            if nm_cell < idx
                                bdry.matid_for_cell{j, i}(nm_cell) = bdry.matid_for_cell{j, i}(idx);
                                bdry.vol_for_cell{j, i}(nm_cell) = bdry.vol_for_cell{j, i}(idx);
                                bdry.mass_for_cell{j, i}(nm_cell) = bdry.mass_for_cell{j, i}(idx);
                                bdry.ener_for_cell{j, i}(nm_cell) = bdry.ener_for_cell{j, i}(idx);
                            end
                            vol_sum = vol_sum + bdry.vol_for_cell{j, i}(nm_cell);
                        end
                    end

                    % Adjust volume to maintain conservation
                    frac = vol_cell / vol_sum;
                    for idx = 1:nm_cell
                        bdry.vol_for_cell{j, i}(idx) = bdry.vol_for_cell{j, i}(idx) * frac;
                    end

                    % Handle original mixed cells
                    if nm_cell > 1
                        nmixcell_int_new = nmixcell_int_new + 1;
                        lsize_mx = lsize_mx + nm_cell;
                    end

                    % Update number of materials in the cell
                    bdry.nmat_for_cell(j, i) = nm_cell;
                end
            end

            bdry.bdry_cell_ragged_2d(mesh, btype_lower, btype_upper);

            nmixcell_new = 0; % new mix cell count
            lsize_mx = 0;     % size of material list in the mixed cells

            for j = 1:ncell_ext(2)
                for i = 1:ncell_ext(1)
                    if bdry.nmat_for_cell(j, i) > 1
                        nmixcell_new = nmixcell_new + 1;
                        lsize_mx = lsize_mx + bdry.nmat_for_cell(j, i);
                    end
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %   mapping velocity

            % allocate 2D double arrays for momentum and velocity
            mom_for_2dnode_new = zeros(nnode_ext(2), nnode_ext(1), dim);  % 2D array for momentum
            mom_for_node = zeros(nnode_ext(2), nnode_ext(1), dim);        % 2D array for momentum

            if if_mom_conserved
                % Initialize the momentum array to zero
                for j = 1:nnode_ext(2)
                    for i = 1:nnode_ext(1)
                        mom_for_node(j, i, :) = 0.0;
                    end
                end

                % Calculate momentum for each node
                for j = 3:nnode_ext(2)-2
                    for i = 3:nnode_ext(1)-2
                        dmass = 0.0;
                        for jc = j-1:j
                            for ic = i-1:i
                                if (mesh.nmat_for_2dcell(jc, ic) == 1) || if_mom_homoginized
                                    dmass = dmass + 0.25 * vol_cell * mesh.rho_for_2dcell(jc, ic);
                                else
                                    coord_node(1) = xl(1) + (i - nbdry) * dx(1);
                                    coord_node(2) = xl(2) + (j - nbdry) * dx(2);

                                    xl_qtr(1) = coord_node(1) + (ic - i) * hdx(1);
                                    xl_qtr(2) = coord_node(2) + (jc - j) * hdx(2);

                                    mixcell = mesh.mixcell_for_2dcell(jc, ic);
                                    [nmat_adv, matid_adv, indice_adv, vol_adv] = mat.intersect2d(xl_qtr, hdx, mixcell);

                                    for idx = 1:nmat_adv
                                        dmass = dmass + vol_adv(idx) * mat.rho_in_mixcell{mixcell}(idx);
                                    end
                                end
                            end
                        end
                        for idx = 1:dim
                            mom_for_node(j, i, idx) = dmass * mesh.vel_for_2dnode(j, i, idx);
                        end
                    end
                end
            else
                % If momentum is not conserved, copy velocities to momentum array
                mom_for_node(:, :, :) = mesh.vel_for_2dnode(:, :, :);
            end

            %   -----------------
            %   |       |       |
            %   |   ....|....   |
            %   |   .   |   .   |
            %   ----.---x---.----
            %   |   .   |   .   |
            %   |   ....|....   |
            %   |       |       |
            %   -----------------

            % Copy the momentum data to the new momentum array
            mom_for_2dnode_new(:, :, :) = mom_for_node(:, :, :);

            if dir == 1
                % Advection in the x-direction
                for j = 3:nnode_ext(2)-2
                    for i = 3:nnode_ext(1)-2
                        % Calculate the distance traveled during advection
                        dist = mesh.vel_for_2dnode(j, i, dir) * dt;
                        frac = abs(dist) / dx(dir);
                        if frac < small, continue; end

                        % Determine the adjacent node in the direction of advection
                        if dist > 0.0
                            ip = i + 1;  % Node to which the momentum is added
                        else
                            ip = i - 1;  % Node to which the momentum is added
                        end

                        % Update momentum at the current and adjacent nodes
                        for idx = 1:dim
                            dmom = frac * mom_for_node(j, i, idx);
                            mom_for_2dnode_new(j, i, idx) = mom_for_2dnode_new(j, i, idx) - dmom;
                            mom_for_2dnode_new(j, ip, idx) = mom_for_2dnode_new(j, ip, idx) + dmom;
                        end
                    end
                end
            elseif dir == 2
                % Advection in the y-direction
                for j = 3:nnode_ext(2)-2
                    for i = 3:nnode_ext(1)-2
                        % Calculate the distance traveled during advection
                        dist = mesh.vel_for_2dnode(j, i, dir) * dt;
                        frac = abs(dist) / dx(dir);
                        if frac < small, continue; end

                        % Determine the adjacent node in the direction of advection
                        if dist > 0.0
                            jp = j + 1;  % Node to which the momentum is added
                        else
                            jp = j - 1;  % Node to which the momentum is added
                        end

                        % Update momentum at the current and adjacent nodes
                        for idx = 1:dim
                            dmom = frac * mom_for_node(j, i, idx);
                            mom_for_2dnode_new(j, i, idx) = mom_for_2dnode_new(j, i, idx) - dmom;
                            mom_for_2dnode_new(jp, i, idx) = mom_for_2dnode_new(jp, i, idx) + dmom;
                        end
                    end
                end
            end

            % put cell variables back to rho_for_2dcell, etc.

            if nmixcell_new > 0
                % Allocate 1D cell arrays for mixed cell data
                ijk_in_mixcell_new = zeros(nmixcell_new, dim);
                nmat_in_mixcell_new = zeros(nmixcell_new, 1);  % Double array for number of materials
                matids_in_mixcell_new = cell(nmixcell_new, 1);
                vf_in_mixcell_new = cell(nmixcell_new, 1);
                rho_in_mixcell_new = cell(nmixcell_new, 1);
                ei_in_mixcell_new = cell(nmixcell_new, 1);
                pres_in_mixcell_new = cell(nmixcell_new, 1);
            else
                ijk_in_mixcell_new = [];
                nmat_in_mixcell_new = [];
                matids_in_mixcell_new = {};
                vf_in_mixcell_new = {};
                rho_in_mixcell_new = {};
                ei_in_mixcell_new = {};
                pres_in_mixcell_new = {};
            end

            % Determine maximum material ID and populate gamma_matid
            matid = max(mat.matids_prob);
            gamma_matid = zeros(matid + 1, 1);
            for i = 1:nmat_prob
                gamma_matid(mat.matids_prob(i)) = gamma_ea_mat(i);
            end

            nmix = 1;
            j_start = [nbdry + 1, 1];
            j_end = [ncell_bdry(2), ncell_ext(2)];
            i_start = [nbdry + 1, 1];
            i_end = [ncell_bdry(1), ncell_ext(1)];

            for pass = 1:2
                jstart = j_start(pass);
                jend = j_end(pass);
                istart = i_start(pass);
                iend = i_end(pass);

                for j = jstart:jend
                    for i = istart:iend
                        if pass == 2 && nbdry < j && j <= ncell_bdry(2) && nbdry < i && i <= ncell_bdry(1)
                            continue;
                        end
                        nm_cell = bdry.nmat_for_cell(j, i);
                        mesh.nmat_for_2dcell(j, i) = nm_cell;
                        mesh.matid_for_2dcell(j, i) = bdry.matid_for_cell{j, i}(1);

                        if nm_cell == 1
                            matid = bdry.matid_for_cell{j, i}(1);
                            mesh.rho_for_2dcell(j, i) = bdry.mass_for_cell{j, i}(1) / vol_cell;
                            mesh.ei_for_2dcell(j, i) = bdry.ener_for_cell{j, i}(1) / vol_cell;
                            if ~is_solid(matid)
                                mesh.pres_for_2dcell(j, i) = (gamma_matid(matid) - 1.0) * mesh.ei_for_2dcell(j, i);
                            else
                                mesh.pres_for_2dcell(j, i) = eos.p_mie_gruneisen(mesh.rho_for_2dcell(j, i), mesh.ei_for_2dcell(j, i));
                            end
                            mesh.pres_for_2dcell(j, i) = max(0.0, mesh.pres_for_2dcell(j, i));
                        else
                            ijk_in_mixcell_new(nmix,1:2) = [i, j];  % Store the indices
                            nmat_in_mixcell_new(nmix) = nm_cell;

                            matids_in_mixcell_new{nmix} = zeros(nm_cell, 1);
                            vf_in_mixcell_new{nmix} = zeros(nm_cell, 1);
                            rho_in_mixcell_new{nmix} = zeros(nm_cell, 1);
                            ei_in_mixcell_new{nmix} = zeros(nm_cell, 1);
                            pres_in_mixcell_new{nmix} = zeros(nm_cell, 1);

                            vol_sum = 0.0;
                            mass_sum = 0.0;
                            ener_sum = 0.0;
                            mesh.pres_for_2dcell(j, i) = 0.0;

                            for idx = 1:nm_cell
                                matid = bdry.matid_for_cell{j, i}(idx);
                                dvol = bdry.vol_for_cell{j, i}(idx);

                                matids_in_mixcell_new{nmix}(idx) = matid;
                                vf_in_mixcell_new{nmix}(idx) = dvol / vol_cell;
                                rho_in_mixcell_new{nmix}(idx) = bdry.mass_for_cell{j, i}(idx) / dvol;
                                ei_in_mixcell_new{nmix}(idx) = bdry.ener_for_cell{j, i}(idx) / dvol;

                                if ~is_solid(matid)
                                    pres_in_mixcell_new{nmix}(idx) = (gamma_matid(matid) - 1.0) * ei_in_mixcell_new{nmix}(idx);
                                else
                                    pres_in_mixcell_new{nmix}(idx) = eos.p_mie_gruneisen(rho_in_mixcell_new{nmix}(idx), ei_in_mixcell_new{nmix}(idx));
                                end
                                pres_in_mixcell_new{nmix}(idx) = max(0.0, pres_in_mixcell_new{nmix}(idx));

                                vol_sum = vol_sum + dvol;
                                mass_sum = mass_sum + bdry.mass_for_cell{j, i}(idx);
                                ener_sum = ener_sum + bdry.ener_for_cell{j, i}(idx);
                                mesh.pres_for_2dcell(j, i) = mesh.pres_for_2dcell(j, i) + (dvol / vol_cell) * pres_in_mixcell_new{nmix}(idx);
                            end
                            mesh.pres_for_2dcell(j, i) = max(0.0, mesh.pres_for_2dcell(j, i));
                            mesh.rho_for_2dcell(j, i) = mass_sum / vol_sum;
                            mesh.ei_for_2dcell(j, i) = ener_sum / vol_sum;
                            nmix = nmix + 1;
                        end
                    end
                end
                if pass == 1
                    assert(nmix - 1 == nmixcell_int_new);
                else
                    assert(nmix - 1 == nmixcell_new);
                end
            end


            %//   apply boundry contion to cell variables
            % mat_set_mix_data(nmixcell_new, nmixcell_int_new, nmat_in_mixcell_new, ijk_in_mixcell_new, matids_in_mixcell_new,
            %                  vf_in_mixcell_new, rho_in_mixcell_new, pres_in_mixcell_new, ei_in_mixcell_new);
            mat.nmixcell = nmixcell_new;
            mat.nmixcell_int = nmixcell_int_new;

            mat.nmat_in_mixcell = nmat_in_mixcell_new;
            mat.ijk_in_mixcell = ijk_in_mixcell_new;
            mat.matids_in_mixcell = matids_in_mixcell_new;
            mat.vf_in_mixcell = vf_in_mixcell_new;
            mat.rho_in_mixcell = rho_in_mixcell_new;
            mat.pres_in_mixcell = pres_in_mixcell_new;
            mat.ei_in_mixcell = ei_in_mixcell_new;

            %   update vel_for_2dnode
            if if_mom_conserved
                % Loop over all nodes within the boundary
                for j = nbdry+1:nnode_bdry(2)
                    for i = nbdry+1:nnode_bdry(1)
                        % Calculate the sum of densities for the surrounding cells
                        rho_sum = 0.25 * sum(sum(mesh.rho_for_2dcell(j-1:j, i-1:i)));

                        % If the sum of densities is too small, set the velocity to zero
                        if rho_sum / rho_average < small
                            mesh.vel_for_2dnode(j, i, 1:dim) = 0.0;
                        else
                            % Update the velocity using the momentum and density sum
                            mesh.vel_for_2dnode(j, i, 1:dim) = mom_for_2dnode_new(j, i, 1:dim) / (vol_cell * rho_sum);
                        end
                    end
                end
            else
                % If momentum is not conserved, copy the momentum values directly to velocity
                mesh.vel_for_2dnode(:, :, :) = mom_for_2dnode_new(:, :, :);
            end

            %   apply boundry contion for vel_for_2dnode
            bdry.bdry_node_2d(mesh, btype_lower, btype_upper);

            % TODO: clear memory

        end

        function [coords_tot, nnode_tot, nface_for_mpoly, nnode_for_face_ea_mpoly, ...
                  nodelist_for_face_ea_mpoly] = build_mpoly3d(obj, mesh, ijk, nm_this_cell)
            % build_mpoly3d Constructs material polyhedra for a 3D mesh.
            %
            % This function constructs material polyhedra within a 3D cell, using volume 
            % fractions and material IDs to define interfaces and faces.
            %
            % Inputs:
            %   obj             - (c_advection) Current c_advection object.
            %   mesh            - (c_mesh) The mesh object.
            %   ijk             - (int array) Indices of the current cell in the 3D grid, size: [3, 1].
            %   nm_this_cell    - (int) Number of materials in the current cell.
            %
            % Outputs:
            %   coords_tot            - (double array) Coordinates of all nodes after reconstruction, 
            %                           size: [3*nnode_tot, 1].
            %   nnode_tot             - (int) Total number of nodes.
            %   nface_for_mpoly       - (int array) Number of faces for each material polygon, size: [nm_this_cell, 1].
            %   nnode_for_face_ea_mpoly - (cell array of int arrays) Number of nodes for each face in 
            %                             each material polygon.
            %   nodelist_for_face_ea_mpoly - (cell array of int arrays) List of node indices for each 
            %                                face in each material polygon.
            %
            % Uses:
            %   mesh.vf_3dmat     - (4D array of double) Volume fractions of materials in 3D cells.
            %   mesh.xl_prob      - (double array) Lower bounds of the problem domain.
            %   mesh.dx           - (double array) Grid spacing in the x, y, z directions.
            %   mat.nmat_prob     - (int) Number of materials in the problem.
            %   mesh.nbdry        - (int) Number of boundary cells.
            %
            % Modifies:
            %   None.
            %
            % Original C declaration:
            % void build_mpoly3d(int nbdry, double *xl_prob, double *dx,
            %                    int nmat, double ****vf_3dmat,
            %                    int *ijk, int nm_this_cell,
            %                    int *nnode_tot, double **coords_tot,
            %                    int *nface_for_mpoly, int ***nnode_for_face_ea_mpoly,
            %                    int ***nodelist_for_face_ea_mpoly)
            global mat;
            
            % Initialize variables
            dim = 3;
            xl = mesh.xl_prob + (ijk - mesh.nbdry - 1) .* mesh.dx;
            dx = mesh.dx;
            nmat = mat.nmat_prob;
            vf_3dmat = mesh.vf_3dmat;
            
            % Temporary buffers for material properties
            nmat_3dsmesh = zeros(3, 3, 3, 'int32');
            matid_3dsmesh = cell(3, 3, 3);
            vf_3dsmesh = cell(3, 3, 3);
        
            % Populate local 3D submesh properties
            for k0 = 0:2
                for j0 = 0:2
                    for i0 = 0:2
                        k = ijk(3) + k0 - 1;
                        j = ijk(2) + j0 - 1;
                        i = ijk(1) + i0 - 1;
                        vf_local = vf_3dmat{k, j, i};
                        matid_local = find(vf_local > 0);
                        nmat_3dsmesh(k0+1, j0+1, i0+1) = numel(matid_local);
                        matid_3dsmesh{k0+1, j0+1, i0+1} = matid_local - 1; % Adjust for MATLAB indexing
                        vf_3dsmesh{k0+1, j0+1, i0+1} = vf_local(matid_local);
                    end
                end
            end
        
            % Call reconstruction function
            [coords_tot, nnode_tot, nodes_for_minterface, nnode_for_face_ea_mpoly, ...
                nodelist_for_face_ea_mpoly] = obj.reconstruct3d_nmat_pagosa(xl, dx, ...
                nmat_3dsmesh, matid_3dsmesh, vf_3dsmesh);
        
        end

        function courant_adv = quantities_crossing_edge_3d(obj, dir, dt, vel_3dface)
            % quantities_crossing_edge_3d Calculates volume, mass, and energy crossing edges in 3D.
            %
            % This method calculates the volume, mass, and energy that cross edges in a
            % 3D grid during advection, accounting for multiple materials and their
            % properties.
            %
            % Inputs:
            %   obj        - (c_advection) The advection object instance.
            %   dir        - (int) The direction of advection (1: x, 2: y, 3: z).
            %   mesh       - (c_mesh) The mesh object containing grid and material data.
            %   dt         - (double) The time step.
            %   vel_3dface - (3D array) Velocity at each face in the 3D grid.
            %
            % Outputs:
            %   courant_adv - (double) Courant factor for advection.
            %
            % Uses:
            %   mesh.vf_3dmat, mesh.rho_3dmat, mesh.ei_3dmat, mesh.ncell, mesh.nbdry,
            %   mesh.dx, mesh.xl_prob, mat.nmat_prob
            %
            % Modifies:
            %   obj.matid_for_edge, obj.dvol_for_edge, obj.dmass_for_edge,
            %   obj.dener_for_edge
            global mesh mat;
            
            % Constants and initialization
            dim = 3;
            ncell = mesh.ncell_prob;
            nbdry = mesh.nbdry_prob;
            dx = mesh.dx;
            xl_prob = mesh.xl_prob;
            nmat = mesh.nmat_mesh;

            vol_cell = prod(dx);

            % Compute extended dimensions and sizes
            ncell_ext = ncell + 2* nbdry;
            ncell_bdry = ncell + nbdry - 1; % cell at the right boundary
            nnode_ext = ncell_ext + 1;
            nnode_bdry = nnode_ext - nbdry - 1; % node at the right boundary
            
            % Define sizes for edges
            sizes_edge = ncell_ext;
            sizes_edge(dir) = nnode_ext(dir);
            lsize_max = prod(sizes_edge) * 2;
        
            obj.matid_for_edge = zeros(lsize_max, 1, 'int32');
            obj.dvol_for_edge = zeros(lsize_max, 1);
            obj.dmass_for_edge = zeros(lsize_max, 1);
            obj.dener_for_edge = zeros(lsize_max, 1);
        
            loc_edge = 1; % Start index for edge quantities
        
            % Compute bounds for loops
            imx = ncell_bdry(1);
            jmx = ncell_bdry(2);
            kmx = ncell_bdry(3);
            if dir == 1, imx = nnode_bdry(1); end
            if dir == 2, jmx = nnode_bdry(2); end
            if dir == 3, kmx = nnode_bdry(3); end
            courant_adv = 0.0;

            for k = nbdry + 1:kmx
                if dir ~= 3
                    xl_cell(3) = xl_prob(3) + (k - nbdry - 1) * dx(3);
                    xl_slab(3) = xl_cell(3);
                    xr_slab(3) = xl_cell(3) + dx(3);
                    inward_norm(3) = 0.0;
                    kl = k;
                else
                    xyz_edge = xl_prob(3) + (k - nbdry - 1) * dx(3);
                    kl = k - 1;
                end
        
                for j = nbdry + 1:jmx
                    if dir ~= 2
                        xl_cell(2) = xl_prob(2) + (j - nbdry - 1) * dx(2);
                        xl_slab(2) = xl_cell(2);
                        xr_slab(2) = xl_cell(2) + dx(2);
                        inward_norm(2) = 0.0;
                        jl = j;
                    else
                        xyz_edge = xl_prob(2) + (j - nbdry - 1) * dx(2);
                        jl = j - 1;
                    end
        
                    for i = nbdry + 1:imx
                        if dir ~= 1
                            xl_cell(1) = xl_prob(1) + (i - nbdry - 1) * dx(1);
                            xl_slab(1) = xl_cell(1);
                            xr_slab(1) = xl_cell(1) + dx(1);
                            inward_norm(1) = 0.0;
                            il = i;
                        else
                            xyz_edge = xl_prob(1) + (i - nbdry - 1) * dx(1);
                            il = i - 1;
                        end
        
                        % Calculate distance and fractional advection
                        dist = vel_3dface(k, j, i) * dt;
                        frac = dist / dx(dir);
                        courant_adv = max(courant_adv, 2.0 * abs(frac));
                        sign_frac = 2*(frac > 0) - 1;
                        
                        % Warn if Courant number exceeds 1
                        if courant_adv >= 1.0
                            warning('Courant number exceeded 1: courant_adv = %e', courant_adv);
                        end
                        
                        % Process advection
                        if abs(frac) > obj.small
                            nm = 0;
                            mats = zeros(nmat, 1);
                            for m = 1:nmat
                                if mesh.vf_3dmat(kl, jl, il, m) >= mat.vfmin
                                    nm = nm + 1;
                                    mats(nm) = m;
                                end
                            end
                            mats = mats(1:nm);
                        
                            if nm == 1
                                % Single material case
                                obj.nmat_for_edge(k, j, i) = 1;
                                dvol = vol_cell * abs(frac);
                                dmass = dvol * mesh.rho_3dmat(kl, jl, il, mats(1));
                                dener = dvol * mesh.ei_3dmat(kl, jl, il, mats(1));
                        
                                obj.matid_for_edge(loc_edge) = mats(1);
                                obj.dvol_for_edge(loc_edge) = dvol;
                                obj.dmass_for_edge(loc_edge) = dmass;
                                obj.dener_for_edge(loc_edge) = dener;
                                loc_edge = loc_edge + 1;
                            else
                                % Multiple material case: build polyhedron and advect
                                if (frac > 0)
                                    xl_slab(dir) = xyz_edge - dist;
                                    xr_slab(dir) = xyz_edge;
                                    xl_cell(dir) = xyz_edge - dx(dir);
                                    inward_norm(dir) = 1.0;
                                    slab_faceid = 2*dir - 1; % 1, 3, 5
                                else
                                    xl_slab(dir) = xyz_edge;
                                    xr_slab(dir) = xyz_edge + abs(dist);
                                    xl_cell(dir) = xyz_edge;
                                    inward_norm(dir) = -1.0;
                                    slab_faceid = 2*dir; % 2, 4, 6
                                end
                        
                                ijk = [il, jl, kl];
                        
                                % Build polyhedron
                                [~, coords_tot, nface_for_mpoly, ...
                                 nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly] = ...
                                 obj.build_mpoly3d(mesh, ijk, nm);
                        
                                % Perform advection
                                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
                                 obj.advect3d(mesh, ijk, xl_cell, nm, mats, ...
                                              xl_slab, xr_slab, inward_norm, slab_faceid, ...
                                              coords_tot, nface_for_mpoly, ...
                                              nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly);
                        
                                % Free temporary arrays
                                delete(coords_tot);
                                delete(nnode_for_face_ea_mpoly);
                                delete(nodelist_for_face_ea_mpoly);
                        
                                % Store results
                                obj.nmat_for_edge(k, j, i) = nmat_adv;
                                obj.matid_for_edge(loc_edge:loc_edge + nmat_adv - 1) = matid_adv(1:nmat_adv);
                                obj.dvol_for_edge(loc_edge:loc_edge + nmat_adv - 1)  = sign_frac*vol_adv(1:nmat_adv);
                                obj.dmass_for_edge(loc_edge:loc_edge + nmat_adv - 1) = sign_frac*mass_adv(1:nmat_adv);
                                obj.dener_for_edge(loc_edge:loc_edge + nmat_adv - 1) = sign_frac*ener_adv(1:nmat_adv);
                                loc_edge = loc_edge + nmat_adv;
                            end
                        end  % if abs(frac) > obj.small

                    end % i
                end % j
            end % k
        end

        function update_materials_after_advection(obj, dir, nmat, sizes_edge)
            % update_materials_after_advection - Update material quantities after advection.
            %
            % This function updates the material volume fractions, densities, and energies
            % in the computational domain after advection is performed in a given direction.
            %
            % Inputs:
            %   obj        - (c_advection) The advection object instance.
            %   dir        - (int) Direction of advection (1, 2, or 3 for x, y, z).
            %   nmat       - (int) Total number of materials in the simulation.
            %   lsize_cell - (int) Total number of cells in the computational domain.
            %   sizes_edge - (array, size [3, 1]) Number of edges in each direction.
            %
            % Uses:
            %   mesh.ncell_prob   - (array, size [3, 1]) Number of cells in each direction.
            %   mesh.dx           - (array, size [3, 1]) Grid spacing in each direction.
            %   mesh.nbdry_prob   - (int) Number of boundary layers around the domain.
            %   mesh.vf_3dmat     - (4D array, size [nz, ny, nx, nmat]) Volume fractions.
            %   mesh.rho_3dmat    - (4D array, size [nz, ny, nx, nmat]) Densities.
            %   mesh.ei_3dmat     - (4D array, size [nz, ny, nx, nmat]) Internal energies.
            %   obj.nmat_for_edge - (3D array, size [nz, ny, nx]) Number of materials at edges.
            %   obj.matid_for_edge - (1D array, dynamic) Material IDs crossing edges.
            %   obj.dvol_for_edge - (1D array, dynamic) Volume changes due to advection.
            %   obj.dmass_for_edge - (1D array, dynamic) Mass changes due to advection.
            %   obj.dener_for_edge - (1D array, dynamic) Energy changes due to advection.
            %
            % Modifies:
            %   obj.vol_for_cell  - (4D array, size [nz, ny, nx, nmat]) Updated cell volumes.
            %   obj.mass_for_cell - (4D array, size [nz, ny, nx, nmat]) Updated cell masses.
            %   obj.ener_for_cell - (4D array, size [nz, ny, nx, nmat]) Updated cell energies.
            %
            % Original C declaration:
            % void update_materials_after_advection(int dir, int nmat, long long lsize_cell, 
            %        int *sizes_edge, int nbdry, int *ncell_ext, int *ncell_bdry, 
            %        int ***nmat_for_edge, double vol_cell, double ****vf_3dmat, 
            %        double ****rho_3dmat, double ****ei_3dmat, double *dvol_for_edge, 
            %        double *dmass_for_edge, double *dener_for_edge, int *matid_for_edge, 
            %        double *****ptr_vol_for_cell, double *****ptr_mass_for_cell, 
            %        double *****ptr_ener_for_cell);

            global mesh mat;
            ncell = mesh.ncell_prob;
            nbdry = mesh.nbdry_prob;
            ncell_ext = ncell + 2*nbdry;  % Add boundary cells
            ncell_bdry = ncell_ext - nbdry - 1;      % Define boundary limits
            loc_for_3dedge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));

            % Loop through all cells in the extended domain
            loc_edge = 1;
            for k = 1:ncell_ext(3)
                for j = 1:ncell_ext(2)
                    for i = 1:ncell_ext(1)
                        % Loop through all materials in the cell
                        for m = 1:nmat
                            % Calculate volume, mass, and energy for each material
                            vfc = vol_cell * mesh.vf_3dmat(k, j, i, m);
                            obj.vol_for_cell(k, j, i, m) = vfc;
                            obj.mass_for_cell(k, j, i, m) = vfc * mesh.rho_3dmat(k, j, i, m);
                            obj.ener_for_cell(k, j, i, m) = vfc * mesh.ei_3dmat(k, j, i, m);
                        end
                        if obj.nmat_for_edge(k, j, i) > 0
                            loc_for_3dedge(k, j, i) = loc_edge;
                            loc_edge = loc_edge + obj.nmat_for_edge(k, j, i);
                        else
                            loc_for_3dedge(k, j, i) = -1;
                        end

                    end
                end
            end

            % Calculate materials after the advection
            for k = 1:sizes_edge(3)
                kc_in = k;
                kc_out = k;
                for j = 1:sizes_edge(2)
                    jc_in = j;
                    jc_out = j;
                    for i = 1:sizes_edge(1)
                        ic_in = i;
                        ic_out = i;
            
                        nm_cross_edge = obj.nmat_for_edge(k, j, i);
                        loc_edge = loc_for_3dedge(k, j, i);
                        if loc_edge < 0
                            continue;
                        end
            
                        % Determine whether materials are entering or leaving the cell
                        if obj.dvol_for_edge(loc_edge) > 0.0
                            if dir == 1  % x-direction
                                ic_out = i - 1;  % Materials to be deducted
                            elseif dir == 2  % y-direction
                                jc_out = j - 1;
                            elseif dir == 3  % z-direction
                                kc_out = k - 1;
                            end
                        else
                            if dir == 1  % x-direction
                                ic_in = i - 1;  % Materials to be added
                            elseif dir == 2  % y-direction
                                jc_in = j - 1;
                            elseif dir == 3  % z-direction
                                kc_in = k - 1;
                            end
                        end
            
                        % Update the incoming materials
                        for idx = 1:nm_cross_edge
                            m = obj.matid_for_edge(loc_edge);
                            obj.vol_for_cell(kc_in, jc_in, ic_in, m) = ...
                            obj.vol_for_cell(kc_in, jc_in, ic_in, m) + abs(obj.dvol_for_edge(loc_edge));
                            obj.mass_for_cell(kc_in, jc_in, ic_in, m) = ...
                            obj.mass_for_cell(kc_in, jc_in, ic_in, m) + abs(obj.dmass_for_edge(loc_edge));
                            obj.ener_for_cell(kc_in, jc_in, ic_in, m) = ...
                            obj.ener_for_cell(kc_in, jc_in, ic_in, m) + abs(obj.dener_for_edge(loc_edge));
                            loc_edge = loc_edge + 1;
                        end
            
                        % Reset loc_edge for outgoing materials
                        loc_edge = loc_for_3dedge(k, j, i);
            
                        % Update the outgoing materials
                        for idx = 1:nm_cross_edge
                            m = obj.matid_for_edge(loc_edge);
                            obj.vol_for_cell(kc_out, jc_out, ic_out, m) = ...
                            obj.vol_for_cell(kc_out, jc_out, ic_out, m) - abs(obj.dvol_for_edge(loc_edge));
                            obj.mass_for_cell(kc_out, jc_out, ic_out, m) = ...
                            obj.mass_for_cell(kc_out, jc_out, ic_out, m) - abs(obj.dmass_for_edge(loc_edge));
                            obj.ener_for_cell(kc_out, jc_out, ic_out, m) = ...
                            obj.ener_for_cell(kc_out, jc_out, ic_out, m) - abs(obj.dener_for_edge(loc_edge));
                            loc_edge = loc_edge + 1;
                        end
                    end
                end
            end

            % Remove small volume fractions
            for k = nbdry:ncell_bdry(3)
                for j = nbdry:ncell_bdry(2)
                    for i = nbdry:ncell_bdry(1)
                        vol_sum = 0.0;
            
                        % Remove small volume fractions below threshold
                        for m = 1:nmat
                            frac = obj.vol_for_cell(k, j, i, m) / vol_cell;
                            if frac < mat.vfmin
                                obj.vol_for_cell(k, j, i, m) = 0.0;
                                obj.mass_for_cell(k, j, i, m) = 0.0;
                                obj.ener_for_cell(k, j, i, m) = 0.0;
                            else
                                vol_sum = vol_sum + obj.vol_for_cell(k, j, i, m);
                            end
                        end
            
                        % Normalize remaining materials to maintain conservation
                        if vol_sum > 0  % Avoid division by zero
                            frac = vol_cell / vol_sum;
                            obj.vol_for_cell(k, j, i, 1:nmat) = ...
                            obj.vol_for_cell(k, j, i, 1:nmat) * frac;
                        end
                    end
                end
            end

        end

        function courant_adv = mapping3d(obj, mesh, mat, dir, dt, ...
                                         btype_lower, btype_upper, ...
                                         is_solid, gamma_ea_mat, ...
                                         vel_face_3d)
            % mapping3d - Perform 3D advection mapping for materials in the mesh.
            %
            % This function handles the advection of materials in a 3D computational
            % mesh. It updates material properties, such as volume fractions,
            % densities, internal energies, and pressures, based on the advection
            % process.
            %
            % Inputs:
            %   obj         - (c_advection) The advection object instance.
            %   mesh        - (c_mesh) The mesh object instance, containing fields
            %                 like ncell_prob, dx, nbdry_prob, and material data.
            %   mat         - (c_mat) Material properties and their corresponding
            %                 methods.
            %   dir         - (int) Direction of advection (1, 2, or 3 for x, y, z).
            %   dt          - (double) Time step for the advection.
            %   btype_lower - (Bdry_Type array) Boundary types for the lower boundary
            %                 (xl, yl, zl).
            %   btype_upper - (Bdry_Type array) Boundary types for the upper boundary
            %                 (xr, yr, zr).
            %   is_solid    - (int array) Solid material indicator array, size: [nmat].
            %   gamma_ea_mat - (double array) Gamma values for each material, size: [nmat].
            %   vel_face_3d - (double array) Velocities at 3D faces, size: [nx, ny, nz].
            %
            % Outputs:
            %   courant_adv - (double) Maximum Courant number observed during the
            %                 advection process.
            %
            % Uses:
            %   mesh.ncell_prob    - (int array) Number of cells in the mesh for each
            %                        dimension, size: [3, 1].
            %   mesh.dx            - (double array) Mesh cell sizes, size: [3, 1].
            %   mesh.nbdry_prob    - (int) Number of boundary layers in the mesh.
            %   mesh.vf_3dmat      - (double 4D array) Volume fractions for materials,
            %                        size: [nz, ny, nx, nmat].
            %   mesh.rho_3dmat     - (double 4D array) Densities for materials, size:
            %                        [nz, ny, nx, nmat].
            %   mesh.ei_3dmat      - (double 4D array) Internal energies for materials,
            %                        size: [nz, ny, nx, nmat].
            %   mesh.pres_3dmat    - (double 4D array) Pressures for materials, size:
            %                        [nz, ny, nx, nmat].
            %   mesh.rho_3dcell    - (double 3D array) Total cell densities, size:
            %                        [nz, ny, nx].
            %   mesh.ei_3dcell     - (double 3D array) Total cell internal energies,
            %                        size: [nz, ny, nx].
            %   mesh.pres_3dcell   - (double 3D array) Total cell pressures, size:
            %                        [nz, ny, nx].
            %   mesh.vel_3dnode    - (double 4D array) Velocities at 3D nodes, size:
            %                        [nz, ny, nx, 3].
            %
            % Modifies:
            %   mesh.vf_3dmat      - (double 4D array) Updated volume fractions for
            %                        materials.
            %   mesh.rho_3dmat     - (double 4D array) Updated densities for materials.
            %   mesh.ei_3dmat      - (double 4D array) Updated internal energies for
            %                        materials.
            %   mesh.pres_3dmat    - (double 4D array) Updated pressures for materials.
            %   mesh.rho_3dcell    - (double 3D array) Updated total cell densities.
            %   mesh.ei_3dcell     - (double 3D array) Updated total cell internal
            %                        energies.
            %   mesh.pres_3dcell   - (double 3D array) Updated total cell pressures.
            %
            % Original C declaration:
            % void mapping3d(int *ncell, int nbdry, double *xl_prob, double *dx, int dir, 
            %                double dt, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
            %                int nmat, int *is_solid, double *gamma_ea_mat,
            %                double ****vf_3dmat, double ****rho_3dmat, double ****ei_3dmat,
            %                double ****pres_3dmat, double ***rho_3dcell, 
            %                double ***ei_3dcell, double ***pres_3dcell, 
            %                double ****vel_3dnode, double ***vel_3dface, 
            %                double *courant_adv);
            global mesh mat bdry;
            
            % Initialize variables
            dim = 3;
            ncell = mesh.ncell_prob;
            dx = mesh.dx;
            nbdry = mesh.nbdry_prob;
        
            % Compute grid dimension variables
            vol_cell = prod(dx);
            ncell_ext = ncell + 2*nbdry;
            ncell_bdry = ncell + nbdry; % cell at the right boundary
            nnode_ext = ncell_ext + 1;
            nnode_bdry = nnode_ext - nbdry; % node at the right boundary
            sizes_edge = ncell_ext;
            sizes_edge(dir) = nnode_ext(dir);
            
     
            % Compute average density
            rho_average = mean(mesh.rho_3dcell(nbdry+1:ncell_bdry(3), ...
                                               nbdry+1:ncell_bdry(2), ...
                                               nbdry+1:ncell_bdry(1)), 'all');
        
            % Allocate and initialize nmat_for_edge
            obj.nmat_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));
            obj.dvol_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));
            obj.dmass_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));
            obj.dener_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));
            obj.matid_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));
        
            % Calculate the volumes, masses, and energies crossing the cell edges
            courant_adv = obj.quantities_crossing_edge_3d(obj, dir, dt, vel_face_3d);
        
            % Update materials after advection
            obj.update_materials_after_advection(dir, nmat, sizes_edge);

            % Mapping velocity

            % Calculate momentum for internal nodes
            mom_for_node = zeros(nnode_ext(3), nnode_ext(2), nnode_ext(1), dim);
            for k = 3:nnode_ext(3)-2
                for j = 3:nnode_ext(2)-2
                    for i = 3:nnode_ext(1)-2
                        dmass = vol_cell * mean(mesh.rho_3dcell(k-1:k,j-1:j,i-1:i),'all');
                        mom_for_node(k, j, i, 1:dim) = dmass * mesh.vel_3dnode(k, j, i, 1:dim);
                    end
                end
            end

            % Copy initial momentum into the new momentum array
            mom_for_3dnode_new = mom_for_node;
            
            % Update momentum based on velocity and advection direction
            for k = 3:nnode_ext(3)-2
                kp = k;
                for j = 3:nnode_ext(2)-2
                    jp = j;
                    for i = 3:nnode_ext(1)-2
                        ip = i;
                        dist = mesh.vel_3dnode(k, j, i, dir) * dt;
                        frac = abs(dist) / dx(dir);
            
                        % Skip if the distance is below the threshold
                        if frac < mesh.small
                            continue;
                        end
            
                        % Determine the node to which the momentum is added
                        if dist > 0.0
                            if dir == 1
                                ip = i + 1; % X-direction
                            elseif dir == 2
                                jp = j + 1; % Y-direction
                            elseif dir == 3
                                kp = k + 1; % Z-direction
                            end
                        else
                            if dir == 1
                                ip = i - 1;
                            elseif dir == 2
                                jp = j - 1;
                            elseif dir == 3
                                kp = k - 1;
                            end
                        end
            
                        % Update momentum for current and neighboring nodes
                        for idx = 1:dim
                            dmom = frac * mom_for_node(k, j, i, idx);
                            mom_for_3dnode_new(k, j, i, idx) = mom_for_3dnode_new(k, j, i, idx) - dmom;
                            mom_for_3dnode_new(kp, jp, ip, idx) = mom_for_3dnode_new(kp, jp, ip, idx) + dmom;
                        end
                    end
                end
            end % k j i

            % Update cell variables for rho, ei, and pres based on the volume fractions
            for k = nbdry+1:ncell_bdry(3)+1
                for j = nbdry+1:ncell_bdry(2)+1
                    for i = nbdry+1:ncell_bdry(1)+1
                        for m = 1:nmat
                            mesh.vf_3dmat(k, j, i, m) = obj.vol_for_cell(k, j, i, m) / mesh.vol_cell;
                            mesh.rho_3dmat(k, j, i, m) = obj.mass_for_cell(k, j, i, m) / ...
                                (mesh.tiny + obj.vol_for_cell(k, j, i, m));
                            mesh.ei_3dmat(k, j, i, m) = obj.ener_for_cell(k, j, i, m) / ...
                                (mesh.tiny + obj.vol_for_cell(k, j, i, m));
                            
                            if is_solid(m)
                                mesh.pres_3dmat(k, j, i, m) = ...
                                    eos.p_mie_gruneisen(mesh.rho_3dmat(k, j, i, m), ...
                                                         mesh.ei_3dmat(k, j, i, m));
                            else 
                                mesh.pres_3dmat(k, j, i, m) = ...
                                    (gamma_ea_mat(m) - 1.0) * mesh.ei_3dmat(k, j, i, m);
                            end
                        end
                    end
                end
            end
            
            % Aggregate material properties into cell-wide properties
            for k = nbdry+1:ncell_bdry(3)
                for j = nbdry+1:ncell_bdry(2)
                    for i = nbdry+1:ncell_bdry(1)
                        vfs = mesh.vf_3dmat(k, j, i, 1:nmat);
                        mesh.rho_3dcell(k, j, i)  = sum(vfs .* mesh.rho_3dmat(k, j, i, 1:nmat));
                        mesh.ei_3dcell(k, j, i)   = sum(vfs .* mesh.ei_3dmat(k, j, i, 1:nmat));
                        mesh.pres_3dcell(k, j, i) = sum(vfs .* mesh.pres_3dmat(k, j, i, 1:nmat));

                    end
                end
            end
            
            % Apply boundary conditions to update the boundary cells
            bdry.bdry_cell_3d(mesh, btype_lower, btype_upper);

            % Update velocity at 3D nodes based on the momentum and density
            for k = nbdry:nnode_bdry(3)
                for j = nbdry:nnode_bdry(2)
                    for i = nbdry:nnode_bdry(1)
                        rho_sum = mean(mesh.rho_3dcell(k-1:k, j-1:j, i-1:i));
                        if rho_sum / rho_average < small
                            mesh.vel_3dnode(k, j, i, 1:dim) = 0.0;
                        else
                            % Update velocity based on momentum and density
                            mesh.vel_3dnode(k, j, i, 1:dim) = ...
                                 mom_for_3dnode_new(k, j, i, 1:dim) / (vol_cell * rho_sum);
                        end
                    end
                end
            end

            % Apply boundary conditions to update the boundary cells
            bdry.bdry_node_3d(mesh, btype_lower, btype_upper);

        end

    end
end
