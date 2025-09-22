classdef SAdvection < handle
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
        if dir == 1
          sizes = [nnode_ext(1), ncell_ext(2), ncell_ext(3)];
        elseif dir == 2
          sizes = [ncell_ext(1), nnode_ext(2), ncell_ext(3)];
        elseif dir == 3
          sizes = [ncell_ext(1), ncell_ext(2), nnode_ext(3)];
        end
        vel_face = zeros(sizes(3), sizes(2), sizes(1));

        if dir == 1
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
        elseif dir == 2
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
        elseif dir == 3
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

    %%
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

      % Get mesh properties
      global tol c_checking h5_checking;
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;

      % Initialize variables
      % ncell_ext = ncell + 2 * nbdry;  % Extended cell size (including boundaries)

      %% Get face velocity (single output)
      vel_face = obj.get_face_velocity(mesh, dir);

      if h5_checking
        if dim == 2
          aux = load_h5debug_2d('debug2d.h5', sprintf("vel_2dface(dir=%d)", dir-1), ncycle);
          assert(norm(vel_face - aux)<tol);
        else
          aux = load_h5debug_3d('debug3d.h5', sprintf("vel_3dface(dir=%d)", dir-1), ncycle);
          assert(norm(vel_face - aux,'fro')<tol);
        end
        printTextWithOK(sprintf('dir=%d - get face velocity',dir),55);
      end

      if dim == 2
        vel_face_2d = vel_face;
        vel_face_3d = [];
      else
        vel_face_2d = [];
        vel_face_3d = vel_face;
      end

      %% mapping
      courant_adv = obj.mapping(mesh, mat, dir, dt, ...
        btype_lower, btype_upper, is_solid, gamma_ea_mat, ...
        vel_face_2d, vel_face_3d,ncycle);

    end

    function courant_adv = mapping(obj, mesh, mat, dir, dt, ...
        btype_lower, btype_upper, is_solid, gamma_ea_mat, ...
        vel_face_2d, vel_face_3d,ncycle)
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
          gamma_ea_mat, vel_face_2d, ncycle);
      elseif dim == 3
        % Call 3D mapping
        courant_adv = obj.mapping3d(dir, dt, ...
          btype_lower, btype_upper, is_solid, ...
          gamma_ea_mat, vel_face_3d, ncycle);
      end
    end

    function courant_adv = mapping2d(obj, mesh, mat, dir, dt, ...
        btype_lower, btype_upper, is_solid, gamma_ea_mat, ...
        vel_face_2d, ncycle)
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

      global bdry small eos tol c_checking h5_checking tiny;
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

      vfmin = mat.vfmin;

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

      %%   calculate the vol, mass, and energy advected across the cell interfaces.
      nmat_for_edge = zeros(sizes_edge(2), sizes_edge(1));  % Initialize the material count at each edge to 0

      % lsize_max = sizes_edge(1) * sizes_edge(2) * 2;  % Maximum size of materials crossing the edges
      % dvol_for_edge  = zeros(lsize_max, 1);  % Array for volumes crossing edges
      % dmass_for_edge = zeros(lsize_max, 1);  % Array for masses crossing edges
      % dener_for_edge = zeros(lsize_max, 1);  % Array for energies crossing edges
      % matid_for_edge = zeros(lsize_max, 1);  % Array for material IDs crossing edges

      loc_edge = 1;  % Initialize current location in dmass_for_edge, etc.
      courant_adv = 0.0; % Initialize courant number and velocity arrays

      %%
      if dir == 1
        % Loop over rows
        for j = nbdry+1:(ncell_bdry(2))
          jc = j;
          xl_cell(2) = xl(2) + (j - nbdry - 1) * dx(2);
          xl_slab(2) = xl_cell(2);  % for the slab advected
          xr_slab(2) = xl_cell(2) + dx(2);
          inward_norm(2) = 0.0;
          nmat = mesh.nmat_mesh;

          for edge = nbdry+1:nnode_bdry(1)
            x_edge = xl(dir) + (edge - nbdry - 1) * dx(dir);

            dist = vel_face_2d(j, edge) * dt;
            frac = dist / dx(dir);
            courant_adv = max(courant_adv, 2.0*abs(frac));
            if courant_adv > 1
              warning("courant_adv = %e in mapping2d", courant_adv);
            end
            clower = edge - 1;
            cupper = edge;

            if j==3 && edge==54 && ncycle==1
              disp("debugging");
            end

            if frac > small
              % [j, edge]
              % Positive advection from clower to cupper
              nm = 0;
              mats = zeros(1, nmat);
              for m = 1:nmat
                if mesh.vf_2dmat(j, clower, m) >= mat.vfmin
                  nm = nm + 1;
                  mats(nm) = m;
                end
              end

              if nm == 1
                % Single material
                m = mats(1);
                % m = find(mats);
                nmat_for_edge(j, edge) = 1;
                dvol = dx(2) * dist;
                dmass = dvol * mesh.rho_2dmat(j, clower, m);
                dener = dvol * mesh.ei_2dmat(j, clower, m);

                matid_for_edge(loc_edge) = m;
                dvol_for_edge(loc_edge) = dvol;
                dmass_for_edge(loc_edge) = dmass;
                dener_for_edge(loc_edge) = dener;
                loc_edge = loc_edge + 1;
              else
                % Define slab boundaries
                xl_slab(dir) = x_edge - abs(dist);
                xr_slab(dir) = x_edge;
                xl_cell(dir) = x_edge - dx(1);
                inward_norm(dir) = 1.0; % Inward norm within the slab

                slab_faceid = 0; % xl-face of slab

                % Define indices for the mixed cell
                ijk = [clower, jc];

                % Build material polygon for the mixed cell
                [nnode_tot, coords_tot, nnode_for_mpoly, nodes_for_mpoly] = ...
                  obj.build_mpoly2d(nbdry, mesh.xl_prob, dx, ...
                  nmat, mesh.vf_2dmat, ijk, nm);

                % Perform advection across the slab

                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
                  mat.advect2d(nmat, ncell, nbdry, mesh.vf_2dmat, mesh.rho_2dmat,...
                  mesh.ei_2dmat, ijk, xl_cell, dx(1:dim), nnode_tot, coords_tot, nm, mats, ...
                  nnode_for_mpoly, nodes_for_mpoly, xl_slab, xr_slab,...
                  inward_norm, slab_faceid, mat.vfmin);
                %%

                % Update material and edge properties
                nmat_for_edge(j, edge) = nmat_adv;
                matid_for_edge(loc_edge + (0:nmat_adv-1)) = matid_adv;
                dvol_for_edge(loc_edge + (0:nmat_adv-1)) = vol_adv;
                dmass_for_edge(loc_edge + (0:nmat_adv-1)) = mass_adv;
                dener_for_edge(loc_edge + (0:nmat_adv-1)) = ener_adv;

                loc_edge = loc_edge + nmat_adv;

              end %nm

            elseif frac < -small
              % [j, j, edge]
              % Negative advection from cupper to clower
              nm = 0;
              mats = zeros(1, nmat);
              for m = 1:nmat
                if vf_2dmat(j, cupper, m) >= vfmin
                  nm = nm + 1;
                  mats(nm) = m;
                end
              end
              adist = abs(dist);
              if nm == 1
                % Single material
                m = mats(1);
                nmat_for_edge(j, edge) = 1;
                dvol = dx(2) * adist;
                dmass = dvol * rho_2dmat(j, cupper, m);
                dener = dvol * ei_2dmat(j, cupper, m);

                loc_edge = loc_edge + 1;
                matid_for_edge(loc_edge) = m;
                dvol_for_edge(loc_edge) = -dvol;
                dmass_for_edge(loc_edge) = -dmass;
                dener_for_edge(loc_edge) = -dener;
              else
                % Define slab boundaries
                xl_slab(dir ) = x_edge;
                xr_slab(dir ) = x_edge + abs(dist);
                xl_cell(dir ) = x_edge;
                inward_norm(dir + 1) = -1.0; % Inward norm within the slab

                slab_faceid = 1; % xr-face of slab

                % Define indices for the mixed cell
                ijk = [cupper, jc];

                % Build material polygon for the mixed cell
                [nnode_tot, coords_tot, nnode_for_mpoly, nodes_for_mpoly] = ...
                  obj.build_mpoly2d(nbdry, mesh.xl_prob, dx, ...
                  nmat, mesh.vf_2dmat, ijk, nm);

                % Perform advection across the slab

                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
                  mat.advect2d(nmat, ncell, nbdry, mesh.vf_2dmat, mesh.rho_2dmat,...
                  mesh.ei_2dmat, ijk, xl_cell, dx(1:dim), nnode_tot, coords_tot, nm, mats, ...
                  nnode_for_mpoly, nodes_for_mpoly, xl_slab, xr_slab,...
                  inward_norm, slab_faceid, mat.vfmin);

                % Update material and edge properties
                nmat_for_edge(j, edge) = nmat_adv;
                matid_for_edge(loc_edge + (1:nmat_adv)) = matid_adv;
                dvol_for_edge(loc_edge + (1:nmat_adv)) = -vol_adv;
                dmass_for_edge(loc_edge + (1:nmat_adv)) = -mass_adv;
                dener_for_edge(loc_edge + (1:nmat_adv)) = -ener_adv;

                loc_edge = loc_edge + nmat_adv;

              end
            end
          end
        end
      elseif dir == 2
        for edge = nbdry:nnode_bdry(2)

          y_edge = mesh.xl_prob(dir) + (edge - nbdry) * dx(dir);

          for i = nbdry:ncell_bdry(1)
            ic = i;
            xl_cell(1) = mesh.xl_prob(1) + (i - nbdry) * dx(1);
            xl_slab(1) = xl_cell(1); % for the slab advected
            xr_slab(1) = xl_cell(1) + dx(1);

            inward_norm(1) = 0.0;
            nmat = mesh.nmat_mesh;

            dist = vel_face_2d(edge, i) * dt;
            frac = dist / dx(dir);
            courant_adv = max(courant_adv, 2.0 * abs(frac));

            if courant_adv >= 1.0
              warning('courant_adv = %e in mapping2d', courant_adv);
            end

            clower = edge - 1;
            cupper = edge;

            if frac > small
              % Positive advection from clower to cupper
              nm = sum(mesh.vf_2dmat(clower, i, :) >= vfmin);
              if nm == 1
                m = find(mesh.vf_2dmat(clower, i, :) >= vfmin, 1);
                dvol = dx(1) * dist;
                dmass = dvol * mesh.rho_2dmat(clower, i, m);
                dener = dvol * mesh.ei_2dmat(clower, i, m);

                matid_for_edge(loc_edge + 1) = m;
                dvol_for_edge(loc_edge + 1) = dvol;
                dmass_for_edge(loc_edge + 1) = dmass;
                dener_for_edge(loc_edge + 1) = dener;
                loc_edge = loc_edge + 1;
              else
                jc = clower;
                xl_slab(2) = y_edge - abs(dist);
                xr_slab(2) = y_edge;
                xl_cell(2) = y_edge - dx(2);
                inward_norm(2) = 1.0; % inward norm within the slab
                slab_faceid = 2; % yl-face of slab

                ijk = [ic, jc];

                % Build material polygon for the mixed cell
                [nnode_tot, coords_tot, nnode_for_mpoly, nodes_for_mpoly] = ...
                  obj.build_mpoly2d(nbdry, mesh.xl_prob, dx, ...
                  nmat, mesh.vf_2dmat, ijk, nm);

                % Perform advection across the slab

                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
                  mat.advect2d(nmat, ncell, nbdry, mesh.vf_2dmat, mesh.rho_2dmat,...
                  mesh.ei_2dmat, ijk, xl_cell, dx(1:dim), nnode_tot, coords_tot, nm, mats, ...
                  nnode_for_mpoly, nodes_for_mpoly, xl_slab, xr_slab,...
                  inward_norm, slab_faceid, mat.vfmin);

                %%

                matid_for_edge(loc_edge + (1:nmat_adv)) = matid_adv;
                dvol_for_edge(loc_edge + (1:nmat_adv)) = vol_adv;
                dmass_for_edge(loc_edge + (1:nmat_adv)) = mass_adv;
                dener_for_edge(loc_edge + (1:nmat_adv)) = ener_adv;
                loc_edge = loc_edge + nmat_adv;
              end
            elseif frac < -small
              % Negative advection from cupper to clower
              nm = sum(mesh.vf_2dmat(cupper, i, :) >= vfmin);
              adist = abs(dist);
              if nm == 1
                m = find(mesh.vf_2dmat(cupper, i, :) >= vfmin, 1);
                dvol = dx(1) * adist;
                dmass = dvol * mesh.rho_2dmat(cupper, i, m);
                dener = dvol * mesh.ei_2dmat(cupper, i, m);

                matid_for_edge(loc_edge + 1) = m;
                dvol_for_edge(loc_edge + 1) = -dvol;
                dmass_for_edge(loc_edge + 1) = -dmass;
                dener_for_edge(loc_edge + 1) = -dener;
                loc_edge = loc_edge + 1;
              else
                jc = cupper;
                xl_slab(2) = y_edge;
                xr_slab(2) = y_edge + abs(dist);
                xl_cell(2) = y_edge;
                inward_norm(2) = -1.0; % inward norm within the slab
                slab_faceid = 3; % yr-face of slab

                ijk = [ic, jc];

                % Build material polygon for the mixed cell
                [nnode_tot, coords_tot, nnode_for_mpoly, nodes_for_mpoly] = ...
                  obj.build_mpoly2d(nbdry, mesh.xl_prob, dx, ...
                  nmat, mesh.vf_2dmat, ijk, nm);

                % Perform advection across the slab

                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
                  mat.advect2d(nmat, ncell, nbdry, mesh.vf_2dmat, mesh.rho_2dmat,...
                  mesh.ei_2dmat, ijk, xl_cell, dx(1:dim), nnode_tot, coords_tot, nm, mats, ...
                  nnode_for_mpoly, nodes_for_mpoly, xl_slab, xr_slab,...
                  inward_norm, slab_faceid, mat.vfmin);

                matid_for_edge(loc_edge + (1:nmat_adv)) = matid_adv;
                dvol_for_edge(loc_edge + (1:nmat_adv)) = -vol_adv;
                dmass_for_edge(loc_edge + (1:nmat_adv)) = -mass_adv;
                dener_for_edge(loc_edge + (1:nmat_adv)) = -ener_adv;
                loc_edge = loc_edge + nmat_adv;
              end
            end % frac
          end %i
        end % edge

      end %dir 2

      %% Allocate 2D array for loc_for_2dedge
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

      %%
      % Initialize 3D arrays to store volume, mass, and energy information
      bdry.vol_for_cell = zeros(ncell_ext(2), ncell_ext(1), nmat);
      bdry.mass_for_cell = zeros(ncell_ext(2), ncell_ext(1), nmat);
      bdry.ener_for_cell = zeros(ncell_ext(2), ncell_ext(1), nmat);

      % Loop through all cells
      for j = 1:ncell_ext(2)
        for i = 1:ncell_ext(1)
          for m = 1:nmat
            % Calculate volumes, masses, and energies
            bdry.vol_for_cell(j, i, m) = vol_cell * mesh.vf_2dmat(j, i, m);
            bdry.mass_for_cell(j, i, m) = bdry.vol_for_cell(j, i, m) * mesh.rho_2dmat(j, i, m);
            bdry.ener_for_cell(j, i, m) = bdry.vol_for_cell(j, i, m) * mesh.ei_2dmat(j, i, m);
          end
        end
      end

      %%  calculate materials after the advection

      if dir == 1
        for j = 1:sizes_edge(2)
          jc = j;
          for i = 1:sizes_edge(1)
            nm_cross_edge = nmat_for_edge(j, i);
            loc_edge = loc_for_2dedge(j, i);
            if loc_edge < 0, continue; end
            % Determine the cells involved in the material exchange
            if dvol_for_edge(loc_edge) > 0.0
              ic_out = i - 1; % The cell where materials are deducted
              ic_in = i;      % The cell where materials are added
            else
              ic_out = i;     % The cell where materials are deducted
              ic_in = i - 1;  % The cell where materials are added
            end

            % Add materials to the incoming cell
            for idx = 1:nm_cross_edge
              m = matid_for_edge(loc_edge);
              bdry.vol_for_cell(jc, ic_in, m) = bdry.vol_for_cell(jc, ic_in, m) + abs(dvol_for_edge(loc_edge));
              bdry.mass_for_cell(jc, ic_in, m) = bdry.mass_for_cell(jc, ic_in, m) + abs(dmass_for_edge(loc_edge));
              bdry.ener_for_cell(jc, ic_in, m) = bdry.ener_for_cell(jc, ic_in, m) + abs(dener_for_edge(loc_edge));
              loc_edge = loc_edge + 1;
            end

            % Deduct materials from the outgoing cell
            loc_edge = loc_for_2dedge(j, i);
            for idx = 1:nm_cross_edge
              m = matid_for_edge(loc_edge);
              bdry.vol_for_cell(jc, ic_out, m) = bdry.vol_for_cell(jc, ic_out, m) - abs(dvol_for_edge(loc_edge));
              bdry.mass_for_cell(jc, ic_out, m) = bdry.mass_for_cell(jc, ic_out, m) - abs(dmass_for_edge(loc_edge));
              bdry.ener_for_cell(jc, ic_out, m) = bdry.ener_for_cell(jc, ic_out, m) - abs(dener_for_edge(loc_edge));
              loc_edge = loc_edge + 1;
            end
          end
        end

      elseif dir == 2
        for j = 1:sizes_edge(2)
          for i = 1:sizes_edge(1)
            ic = i;
            nm_cross_edge = nmat_for_edge(j, i);
            loc_edge = loc_for_2dedge(j, i);
            if loc_edge < 0, continue; end

            % Determine the cells involved in the material exchange
            if dvol_for_edge(loc_edge) > 0.0
              jc_out = j - 1; % The cell where materials are deducted
              jc_in = j;      % The cell where materials are added
            else
              jc_out = j;     % The cell where materials are deducted
              jc_in = j - 1;  % The cell where materials are added
            end

            % Add materials to the incoming cell
            for idx = 1:nm_cross_edge
              m = matid_for_edge(loc_edge);
              bdry.vol_for_cell(jc_in, ic, m) = bdry.vol_for_cell(jc_in, ic, m) + abs(dvol_for_edge(loc_edge));
              bdry.mass_for_cell(jc_in, ic, m) = bdry.mass_for_cell(jc_in, ic, m) + abs(dmass_for_edge(loc_edge));
              bdry.ener_for_cell(jc_in, ic, m) = bdry.ener_for_cell(jc_in, ic, m) + abs(dener_for_edge(loc_edge));
              loc_edge = loc_edge + 1;
            end

            % Deduct materials from the outgoing cell
            loc_edge = loc_for_2dedge(j, i);
            for idx = 1:nm_cross_edge
              m = matid_for_edge(loc_edge);
              bdry.vol_for_cell(jc_out, ic, m) = bdry.vol_for_cell(jc_out, ic, m) - abs(dvol_for_edge(loc_edge));
              bdry.mass_for_cell(jc_out, ic, m) = bdry.mass_for_cell(jc_out, ic, m) - abs(dmass_for_edge(loc_edge));
              bdry.ener_for_cell(jc_out, ic, m) = bdry.ener_for_cell(jc_out, ic, m) - abs(dener_for_edge(loc_edge));
              loc_edge = loc_edge + 1;
            end
          end
        end
      end

      %%  take out small volume fraction
      for j = nbdry+1:ncell_bdry(2)
        for i = nbdry+1:ncell_bdry(1)
          vol_sum = 0.0;
          for m = 1:nmat
            frac = bdry.vol_for_cell(j, i, m) / vol_cell;
            if frac < vfmin
              bdry.vol_for_cell(j, i, m) = 0.0;
              bdry.mass_for_cell(j, i, m) = 0.0;
              bdry.ener_for_cell(j, i, m) = 0.0;
            else
              vol_sum = vol_sum + bdry.vol_for_cell(j, i, m);
            end
          end
          frac = vol_cell / vol_sum;
          for m = 1:nmat
            bdry.vol_for_cell(j, i, m) = bdry.vol_for_cell(j, i, m) * frac;
            % Uncomment the following lines if you need to scale mass and energy:
            % bdry.mass_for_cell(j, i, m) = bdry.mass_for_cell(j, i, m) * frac;
            % bdry.ener_for_cell(j, i, m) = bdry.ener_for_cell(j, i, m) * frac;
          end
        end
      end

      %% %%%%%%%%%%%%%%%%%%%%%%%%
      %   mapping velocity

      % allocate 2D double arrays for momentum and velocity
      % Initialize 3D arrays for momentum and velocity
      % mom_for_2dnode_new = zeros(nnode_ext(2), nnode_ext(1), dim);  % New momentum array
      mom_for_node = zeros(nnode_ext(2), nnode_ext(1), dim);        % Current momentum array

      % Zero initialize momentum array
      mom_for_node(:, :, :) = 0.0;

      % Loop over nodes within the extended boundary
      for j = 3:nnode_ext(2)-2
        for i = 3:nnode_ext(1)-2
          dmass = 0.0;

          % Compute the sum of densities from surrounding cells
          for jc = j-1:j
            for ic = i-1:i
              dmass = dmass + (0.25 * vol_cell * mesh.rho_for_2dcell(jc, ic));
            end
          end

          % Update momentum for each dimension
          for idx = 1:dim
            mom_for_node(j, i, idx) = dmass * mesh.vel_for_2dnode(j, i, idx);
          end
        end
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

      %% Copy momentum data to the new momentum array
      mom_for_2dnode_new = mom_for_node;

      if dir == 1
        % Advection in the x-direction
        for j = 3:nnode_ext(2)-2
          for i = 3:nnode_ext(1)-2
            dist = mesh.vel_for_2dnode(j, i, dir) * dt;
            frac = abs(dist) / dx(dir);
            if frac < small
              continue;
            end

            % Determine the node to which momentum will be added
            if dist > 0.0
              ip = i + 1; % Rightward advection
            else
              ip = i - 1; % Leftward advection
            end

            % Update momentum for each dimension
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
            dist = mesh.vel_for_2dnode(j, i, dir) * dt;
            frac = abs(dist) / dx(dir);
            if frac < small
              continue;
            end

            % Determine the node to which momentum will be added
            if dist > 0.0
              jp = j + 1; % Upward advection
            else
              jp = j - 1; % Downward advection
            end

            % Update momentum for each dimension
            for idx = 1:dim
              dmom = frac * mom_for_node(j, i, idx);
              mom_for_2dnode_new(j, i, idx) = mom_for_2dnode_new(j, i, idx) - dmom;
              mom_for_2dnode_new(jp, i, idx) = mom_for_2dnode_new(jp, i, idx) + dmom;
            end
          end
        end
      end

      %% put cell variables back to rho_for_2dcell, etc.
      % Update per-material properties in each cell

      for j = nbdry+1:ncell_bdry(2)
        for i = nbdry+1:ncell_bdry(1)
          for m = 1:nmat
            % Update volume fraction, density, and internal energy
            mesh.vf_2dmat(j, i, m) = bdry.vol_for_cell(j, i, m) / vol_cell;
            mesh.rho_2dmat(j, i, m) = bdry.mass_for_cell(j, i, m) / (tiny +bdry.vol_for_cell(j, i, m));
            mesh.ei_2dmat(j, i, m) = bdry.ener_for_cell(j, i, m) / (tiny + bdry.vol_for_cell(j, i, m));

            % Update pressure
            if ~is_solid(m)
              mesh.pres_2dmat(j, i, m) = (gamma_ea_mat(m) - 1.0) * mesh.ei_2dmat(j, i, m);
            else
              mesh.pres_2dmat(j, i, m) = eos.p_mie_gruneisen(mesh.rho_2dmat(j, i, m), mesh.ei_2dmat(j, i, m));
            end
          end
        end
      end

      % Aggregate cell properties from materials
      for j = nbdry+1:ncell_bdry(2)
        for i = nbdry+1:ncell_bdry(1)
          % Initialize cell properties
          mesh.rho_for_2dcell(j, i) = 0.0;
          mesh.ei_for_2dcell(j, i) = 0.0;
          mesh.pres_for_2dcell(j, i) = 0.0;

          % Aggregate contributions from materials
          for m = 1:nmat
            mesh.rho_for_2dcell(j, i) = mesh.rho_for_2dcell(j, i) + mesh.vf_2dmat(j, i, m) * mesh.rho_2dmat(j, i, m);
            mesh.ei_for_2dcell(j, i) = mesh.ei_for_2dcell(j, i) + mesh.vf_2dmat(j, i, m) * mesh.ei_2dmat(j, i, m);
            mesh.pres_for_2dcell(j, i) = mesh.pres_for_2dcell(j, i) + mesh.vf_2dmat(j, i, m) * mesh.pres_2dmat(j, i, m);
          end
        end
      end

      %% apply boundary condition
      bdry.bdry_cell_2d(mesh,btype_lower, btype_upper);

      %% update vel_for_2dnde
      for j = nbdry+1:nnode_bdry(2)
        for i = nbdry+1:nnode_bdry(1)
          rho_sum = 0.0;

          % Sum the densities of surrounding cells
          for jc = j-1:j
            for ic = i-1:i
              rho_sum = rho_sum + 0.25 * mesh.rho_for_2dcell(jc, ic);
            end
          end

          % Check if the density ratio is below the threshold
          if rho_sum / rho_average < small
            % Set velocities to zero if density ratio is too small
            mesh.vel_for_2dnode(j, i, 1:dim) = 0.0;
          else
            % Update velocity values
            for idx = 1:dim
              mesh.vel_for_2dnode(j, i, idx) = mom_for_2dnode_new(j, i, idx) / (vol_cell * rho_sum);
            end
          end
        end
      end

      bdry.bdry_node_2d(mesh,btype_lower, btype_upper);

    end %mapping2d

    function courant_adv = mapping3d(obj, dir, dt, ...
        btype_lower, btype_upper, ...
        is_solid, gamma_ea_mat, ...
        vel_face_3d, ncycle)
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
      %   ncycle      - (int) Current time iteration cycle.
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
      %   mesh.rho_for_3dcell - (double 3D array) Total cell densities, size:
      %                        [nz, ny, nx].
      %   mesh.ei_for_3dcell - (double 3D array) Total cell internal energies,
      %                        size: [nz, ny, nx].
      %   mesh.pres_for_3dcell - (double 3D array) Total cell pressures, size:
      %                        [nz, ny, nx].
      %   mesh.vel_for_3dnode    - (double 4D array) Velocities at 3D nodes, size:
      %                        [nz, ny, nx, 3].
      %
      % Modifies:
      %   mesh.vf_3dmat      - (double 4D array) Updated volume fractions for
      %                        materials.
      %   mesh.rho_3dmat     - (double 4D array) Updated densities for materials.
      %   mesh.ei_3dmat      - (double 4D array) Updated internal energies for
      %                        materials.
      %   mesh.pres_3dmat    - (double 4D array) Updated pressures for materials.
      %   mesh.rho_for_3dcell- (double 3D array) Updated total cell densities.
      %   mesh.ei_for_3dcell - (double 3D array) Updated total cell internal
      %                        energies.
      %   mesh.pres_for_3dcell- (double 3D array) Updated total cell pressures.
      %
      % Original C declaration:
      % void mapping3d(int *ncell, int nbdry, double *xl_prob, double *dx, int dir,
      %                double dt, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                int nmat, int *is_solid, double *gamma_ea_mat,
      %                double ****vf_3dmat, double ****rho_3dmat, double ****ei_3dmat,
      %                double ****pres_3dmat, double ***rho_for_3dcell,
      %                double ***ei_for_3dcell, double ***pres_for_3dcell,
      %                double ****vel_for_3dnode, double ***vel_3dface,
      %                double *courant_adv);
      global mesh bdry tiny h5_checking tol;

      % Initialize variables
      dim = 3;
      ncell = mesh.ncell_prob;
      dx = mesh.dx;
      nbdry = mesh.nbdry_prob;
      nmat = mesh.nmat_mesh;

      % Compute grid dimension variables
      vol_cell = prod(dx);
      ncell_ext = ncell + 2*nbdry;
      ncell_bdry = ncell + nbdry; % cell at the right boundary
      nnode_ext = ncell_ext + 1;
      nnode_bdry = nnode_ext - nbdry; % node at the right boundary
      sizes_edge = ncell_ext;
      sizes_edge(dir) = nnode_ext(dir);


      %% Compute average density
      rho_average = mean(mesh.rho_for_3dcell(nbdry+1:ncell_bdry(3), ...
        nbdry+1:ncell_bdry(2), ...
        nbdry+1:ncell_bdry(1)), 'all');

      %
      %% Allocate and initialize nmat_for_edge
      obj.nmat_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));
      obj.dvol_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));
      obj.dmass_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));
      obj.dener_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));
      obj.matid_for_edge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));

      % Calculate the volumes, masses, and energies crossing the cell edges
      courant_adv = obj.quantities_crossing_edge_3d(dir, dt, vel_face_3d);

      if h5_checking
        aux = load_h5debug_3d('debug3d.h5',sprintf("nmat_for_edge(dir=%d)",(dir-1)),ncycle);
        assert(norm(obj.nmat_for_edge - double(aux),'fro')<tol);
        printTextWithOK(sprintf('dir=%d - nmat_for_edge ',dir),55);

        aux = load_h5debug_3d('debug3d.h5',sprintf("dvol_for_edge(dir=%d)",(dir-1)),ncycle);
        assert(norm(squeeze(obj.dvol_for_edge(:)) - double(squeeze(aux(1:numel(obj.dvol_for_edge)))),'fro')<tol);
        printTextWithOK(sprintf('dir=%d - dvol_for_edge ',dir),55);

        aux = load_h5debug_3d('debug3d.h5',sprintf("dmass_for_edge(dir=%d)",(dir-1)),ncycle);
        assert(norm(obj.dmass_for_edge(:) - double(squeeze(aux(1:numel(obj.dmass_for_edge)))),'fro')<tol);
        printTextWithOK(sprintf('dir=%d - dmass_for_edge ',dir),55);

        aux = load_h5debug_3d('debug3d.h5',sprintf("dener_for_edge(dir=%d)",(dir-1)),ncycle);
        assert(norm(obj.dener_for_edge(:) - double(squeeze(aux(1:numel(obj.dener_for_edge)))),'fro')<tol);
        printTextWithOK(sprintf('dir=%d - dener_for_edge ',dir),55);
      end

      %% Update materials after advection
      obj.update_materials_after_advection(dir, nmat, vol_cell, sizes_edge, ncycle);
      

      %% Mapping velocity
      % Calculate momentum for internal nodes
      mom_for_node = zeros(nnode_ext(3), nnode_ext(2), nnode_ext(1), dim);
      for k = 2:nnode_ext(3)-1
        for j = 2:nnode_ext(2)-1
          for i = 2:nnode_ext(1)-1
            dmass = vol_cell * mean(mesh.rho_for_3dcell(k-1:k,j-1:j,i-1:i),'all');
            mom_for_node(k, j, i, 1:dim) = dmass * mesh.vel_for_3dnode(k, j, i, 1:dim);
          end
        end
      end
      
      if h5_checking
        aux = load_h5debug_4d('debug3d.h5',sprintf("mom_for_node(dir=%d)",(dir-1)),ncycle);
        assert(norm(mom_for_node - aux,'fro')<tol);
        printTextWithOK(sprintf('dir=%d - mom_for_node ',dir),55);
      end

      %% Copy initial momentum into the new momentum array
      mom_for_3dnode_new = mom_for_node;

      % Update momentum based on velocity and advection direction
      % for k = 3:nnode_ext(3)-2
      %   kp = k;
      %   for j = 3:nnode_ext(2)-2
      %     jp = j;
      %     for i = 3:nnode_ext(1)-2
      %       ip = i;
      %       dist = mesh.vel_for_3dnode(k, j, i, dir) * dt;
      %       frac = abs(dist) / dx(dir);
      % 
      %       % Skip if the distance is below the threshold
      %       if frac < obj.small
      %         continue;
      %       end
      % 
      %       % Determine the node to which the momentum is added
      %       if dist > 0.0
      %         if dir == 1
      %           ip = i + 1; % X-direction
      %         elseif dir == 2
      %           jp = j + 1; % Y-direction
      %         elseif dir == 3
      %           kp = k + 1; % Z-direction
      %         end
      %       else
      %         if dir == 1
      %           ip = i - 1;
      %         elseif dir == 2
      %           jp = j - 1;
      %         elseif dir == 3
      %           kp = k - 1;
      %         end
      %       end
      % 
      %       % Update momentum for current and neighboring nodes
      %       for idx = 1:dim
      %         dmom = frac * mom_for_node(k, j, i, idx);
      %         mom_for_3dnode_new(k, j, i, idx) = mom_for_3dnode_new(k, j, i, idx) - dmom;
      %         mom_for_3dnode_new(kp, jp, ip, idx) = mom_for_3dnode_new(kp, jp, ip, idx) + dmom;
      %       end
      %     end
      %   end
      % end % k j i
      
      % New code
      for k = 2:nnode_ext(3)-1
        for j = 2:nnode_ext(2)-1
          for i = 2:nnode_ext(1)-1
            dist = mesh.vel_for_3dnode(k, j, i, dir) * dt;
            frac = abs(dist) / dx(dir);

            % Skip if below threshold
            if frac < obj.small
              continue;
            end

            % Initialize ip/jp/kp to current node
            ip = i;
            jp = j;
            kp = k;

            % Determine target node (positive or negative direction)
            if dist > 0.0
              if dir == 1
                ip = i + 1;  % X-direction
              elseif dir == 2
                jp = j + 1;  % Y-direction
              elseif dir == 3
                kp = k + 1;  % Z-direction
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

            % Perform momentum update
            for idx = 1:dim
              dmom = frac * mom_for_node(k, j, i, idx);
              mom_for_3dnode_new(k, j, i, idx)  = mom_for_3dnode_new(k, j, i, idx)  - dmom;
              mom_for_3dnode_new(kp, jp, ip, idx) = mom_for_3dnode_new(kp, jp, ip, idx) + dmom;
            end
          end
        end
      end
      %% check mom_for_3dnode_new
      if h5_checking
        aux = load_h5debug_4d('debug3d.h5',...
          sprintf("mom_for_3dnode_new(dir=%d)",(dir-1)),ncycle);
        assert(norm(mom_for_3dnode_new - aux,'fro')<tol);
        printTextWithOK(sprintf('dir=%d - mom_for_node_new',dir),55);
        aux = load_h5debug_4d('debug3d.h5', ...
          sprintf("before_put_cell:vf_3dmat(dir=%d)", dir-1), ncycle);
        assert(norm(mesh.vf_3dmat - aux, 'fro') < tol);
      end

      %% Update cell variables for rho, ei, and pres based on the volume fractions
      for k = nbdry+1:ncell_bdry(3)
        for j = nbdry+1:ncell_bdry(2)
          for i = nbdry+1:ncell_bdry(1)
            for m = 1:nmat
              % mesh.vf_3dmat(k, j, i, m) = min(1.0, obj.vol_for_cell(k, j, i, m) / vol_cell);
              mesh.vf_3dmat(k, j, i, m) = obj.vol_for_cell(k, j, i, m) / vol_cell;
              mesh.rho_3dmat(k, j, i, m) = obj.mass_for_cell(k, j, i, m) / ...
                (tiny + obj.vol_for_cell(k, j, i, m));
              mesh.ei_3dmat(k, j, i, m) = obj.ener_for_cell(k, j, i, m) / ...
                (tiny + obj.vol_for_cell(k, j, i, m));

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
      
      %% check put cell back -- mat variables
      if h5_checking
        % dir is assumed to be in MATLAB indexing (1-based)
        dir_c = dir - 1;

        % === 4D variables ===

        aux = load_h5debug_4d('debug3d.h5', sprintf("put_cell_back:pres_3dmat(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.pres_3dmat - aux, 'fro') < tol);
        % printTextWithOK(sprintf('dir=%d - pres_3dmat', dir), 55);

        aux = load_h5debug_4d('debug3d.h5', sprintf("put_cell_back:ei_3dmat(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.ei_3dmat - aux, 'fro') < tol);
        % printTextWithOK(sprintf('dir=%d - ei_3dmat', dir), 55);

        aux = load_h5debug_4d('debug3d.h5', sprintf("put_cell_back:rho_3dmat(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.rho_3dmat - aux, 'fro') < tol);
        % printTextWithOK(sprintf('dir=%d - rho_3dmat', dir), 55);

        aux = load_h5debug_4d('debug3d.h5', sprintf("put_cell_back:vf_3dmat(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.vf_3dmat - aux, 'fro') < tol);
        % printTextWithOK(sprintf('dir=%d - vf_3dmat', dir), 55);
      end

      %% Aggregate material properties into cell-wide properties
      for k = nbdry+1:ncell_bdry(3)
        for j = nbdry+1:ncell_bdry(2)
          for i = nbdry+1:ncell_bdry(1)
            vfs = mesh.vf_3dmat(k, j, i, 1:nmat);
            mesh.rho_for_3dcell(k, j, i)  = sum(vfs .* mesh.rho_3dmat(k, j, i, 1:nmat));
            mesh.ei_for_3dcell(k, j, i)   = sum(vfs .* mesh.ei_3dmat(k, j, i, 1:nmat));
            mesh.pres_for_3dcell(k, j, i) = sum(vfs .* mesh.pres_3dmat(k, j, i, 1:nmat));

          end
        end
      end
      
      %% check put cell back section
      if h5_checking
        % dir is assumed to be in MATLAB indexing (1-based)
        dir_c = dir - 1;
        % === 3D  cell variables ===

        aux = load_h5debug_3d('debug3d.h5', sprintf("put_cell_back:pres_3dcell(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.pres_for_3dcell - aux, 'fro') < tol);
        % printTextWithOK(sprintf('dir=%d - pres_3dcell', dir), 55);

        aux = load_h5debug_3d('debug3d.h5', sprintf("put_cell_back:ei_3dcell(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.ei_for_3dcell - aux, 'fro') < tol);
        % printTextWithOK(sprintf('dir=%d - ei_3dcell', dir), 55);

        aux = load_h5debug_3d('debug3d.h5', sprintf("put_cell_back:rho_3dcell(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.rho_for_3dcell - aux, 'fro') < tol);

        printTextWithOK(sprintf('dir=%d - put cell back step', dir), 55);
      end



      %% Apply boundary conditions to update the boundary cells
      bdry.bdry_cell_3d(mesh, btype_lower, btype_upper);
      
      %%
      if h5_checking
        % dir_c = dir - 1;  % convert MATLAB dir (1-based) to C dir (0-based)
        % === 4D variables ===
        aux = load_h5debug_4d('debug3d.h5', sprintf("end_of_advec:pres_3dmat(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.pres_3dmat - aux, 'fro') < tol);
        printTextWithOK(sprintf('dir=%d - end_of_advec:pres_3dmat', dir), 55);

        aux = load_h5debug_4d('debug3d.h5', sprintf("end_of_advec:ei_3dmat(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.ei_3dmat - aux, 'fro') < tol);
        printTextWithOK(sprintf('dir=%d - end_of_advec:ei_3dmat', dir), 55);

        aux = load_h5debug_4d('debug3d.h5', sprintf("end_of_advec:rho_3dmat(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.rho_3dmat - aux, 'fro') < tol);
        printTextWithOK(sprintf('dir=%d - end_of_advec:rho_3dmat', dir), 55);

        aux = load_h5debug_4d('debug3d.h5', sprintf("end_of_advec:vf_3dmat(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.vf_3dmat - aux, 'fro') < tol);
        printTextWithOK(sprintf('dir=%d - end_of_advec:vf_3dmat', dir), 55);

        % === 3D variables ===

        aux = load_h5debug_3d('debug3d.h5', sprintf("end_of_advec:pres_3dcell(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.pres_for_3dcell - aux, 'fro') < tol);
        printTextWithOK(sprintf('dir=%d - end_of_advec:pres_3dcell', dir), 55);

        aux = load_h5debug_3d('debug3d.h5', sprintf("end_of_advec:ei_3dcell(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.ei_for_3dcell - aux, 'fro') < tol);
        printTextWithOK(sprintf('dir=%d - end_of_advec:ei_3dcell', dir), 55);

        aux = load_h5debug_3d('debug3d.h5', sprintf("end_of_advec:rho_3dcell(dir=%d)", dir_c), ncycle);
        assert(norm(mesh.rho_for_3dcell - aux, 'fro') < tol);
        printTextWithOK(sprintf('dir=%d - end_of_advec:rho_3dcell', dir), 55);
      end
      %% Update velocity at 3D nodes based on the momentum and density
      for k = nbdry:nnode_bdry(3)
        for j = nbdry:nnode_bdry(2)
          for i = nbdry:nnode_bdry(1)
            rho_sum = mean(mesh.rho_for_3dcell(k-1:k, j-1:j, i-1:i), 'all');
            if rho_sum / rho_average < obj.small
              mesh.vel_for_3dnode(k, j, i, 1:dim) = 0.0;
            else
              % Update velocity based on momentum and density
              mesh.vel_for_3dnode(k, j, i, 1:dim) = ...
                mom_for_3dnode_new(k, j, i, 1:dim) / (vol_cell * rho_sum);
            end
          end
        end
      end

      %% Apply boundary conditions to update the boundary cells
      bdry.bdry_node_3d(mesh, btype_lower, btype_upper);
      
      %% check node variables
      if h5_checking
        aux = load_h5debug_4d('debug3d.h5', ...
          sprintf("end_of_advec:vel_3dnode(dir=%d)", dir-1), ncycle);
        assert(norm(mesh.vel_for_3dnode - aux, 'fro') < tol);
        printTextWithOK(sprintf('dir=%d - end_of_advec:rho_3dcell', dir), 55);
      end
    end
    
    %%
    function [nnode_tot, coords_tot, nface_for_mpoly, nnode_for_face_ea_mpoly, ...
        nodelist_for_face_ea_mpoly] = build_mpoly3d(~, mesh, ijk)
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
      global vof3d;

      % Initialize variables
      dim = 3;
      xl = mesh.xl_prob + (ijk - mesh.nbdry_prob - 1) .* mesh.dx;
      dx = mesh.dx;
      nmat = mesh.nmat_mesh;

      % Temporary buffers for material properties
      nmat_3dsmesh = zeros(3, 3, 3, 'int32');
      matid_3dsmesh = -ones(3, 3, 3, nmat);
      vf_3dsmesh = zeros(3, 3, 3, nmat);

      % Populate local 3D submesh properties
      for k0 = 1:3
        k = ijk(3) + k0 - 2;
        for j0 = 1:3
          j = ijk(2) + j0 - 2;
          for i0 = 1:3
            i = ijk(1) + i0 - 2;
            vf_local = mesh.vf_3dmat(k, j, i, 1:nmat);
            nm = 0;
            for m = 1:nmat
              if vf_local(m) > 0
                nm = nm + 1;
                matid_3dsmesh(k0, j0, i0, nm) = m;
                vf_3dsmesh(k0, j0, i0, nm) = vf_local(m);
              end
            end
            nmat_3dsmesh(k0, j0, i0) = nm;
          end
        end
      end

      % Call reconstruction function
      [nnode_tot, coords_tot, ~, ~, nface_for_mpoly, nnode_for_face_ea_mpoly, ...
        nodelist_for_face_ea_mpoly] = vof3d.reconstruct3d_nmat_pagosa(xl, dx, ...
        nmat_3dsmesh, matid_3dsmesh, vf_3dsmesh);

    end
    
    %%
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

      global h5_checking

      %% Constants and initialization
      dim = 3;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      dx = mesh.dx;
      xl_prob = mesh.xl_prob;
      nmat = mesh.nmat_mesh;

      vol_cell = prod(dx);

      % Compute extended dimensions and sizes
      ncell_ext = ncell + 2* nbdry;
      ncell_bdry = ncell + nbdry; % cell at the right boundary
      nnode_ext = ncell_ext + 1;
      nnode_bdry = nnode_ext - nbdry; % node at the right boundary

      % Define sizes for edges
      sizes_edge = ncell_ext;
      sizes_edge(dir) = nnode_ext(dir);
      lsize_max = prod(sizes_edge) * 2;
      %
      obj.matid_for_edge = zeros(lsize_max, 1, 'int32');
      obj.dvol_for_edge = zeros(lsize_max, 1);
      obj.dmass_for_edge = zeros(lsize_max, 1);
      obj.dener_for_edge = zeros(lsize_max, 1);

      if h5_checking
        obj.matid_for_edge = [];
        obj.dvol_for_edge = [];
        obj.dmass_for_edge = [];
        obj.dener_for_edge = [];
      end


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

            %%
            % if h5_checking
            %   if loc_edge == 7
            %     display(loc_edge);
            %     keyboard
            %   end
            % end
            %% Process advection
            if frac > obj.small
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
                dvol = sign_frac*vol_cell * abs(frac);
                dmass = dvol * mesh.rho_3dmat(kl, jl, il, mats(1));
                dener = dvol * mesh.ei_3dmat(kl, jl, il, mats(1));

                obj.matid_for_edge(loc_edge) = mats(1);
                obj.dvol_for_edge(loc_edge) = dvol;
                obj.dmass_for_edge(loc_edge) = dmass;
                obj.dener_for_edge(loc_edge) = dener;
                loc_edge = loc_edge + 1;
              else
                % Multiple material case: build polyhedron and advect
                xl_slab(dir) = xyz_edge - dist;
                xr_slab(dir) = xyz_edge;
                xl_cell(dir) = xyz_edge - dx(dir);
                inward_norm(dir) = 1.0;
                slab_faceid = 2*dir - 1; % 1, 3, 5

                ijk = [il, jl, kl];

                %% Build polyhedron
                [nnode_tot, coords_tot, nface_for_mpoly, ...
                  nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly] = ...
                  obj.build_mpoly3d(mesh, ijk);

                % Perform advection
                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
                  mat.advect3d(mesh, ijk, xl_cell, dx, nm, mats, ...
                  xl_slab, xr_slab, inward_norm, slab_faceid, ...
                  nnode_tot, coords_tot, nface_for_mpoly, ...
                  nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly);

                % Free temporary arrays
                coords_tot = [];
                nnode_for_face_ea_mpoly = [];
                nodelist_for_face_ea_mpoly = [];

                % Store results
                obj.nmat_for_edge(k, j, i) = nmat_adv;
                obj.matid_for_edge(loc_edge:loc_edge + nmat_adv - 1) = matid_adv(1:nmat_adv);
                obj.dvol_for_edge(loc_edge:loc_edge + nmat_adv - 1)  = sign_frac*vol_adv(1:nmat_adv);
                obj.dmass_for_edge(loc_edge:loc_edge + nmat_adv - 1) = sign_frac*mass_adv(1:nmat_adv);
                obj.dener_for_edge(loc_edge:loc_edge + nmat_adv - 1) = sign_frac*ener_adv(1:nmat_adv);
                loc_edge = loc_edge + nmat_adv;

              end
            elseif frac<-obj.small
              nm = 0;
              mats = zeros(nmat, 1);
              for m = 1:nmat
                if mesh.vf_3dmat(k, j, i, m) >= mat.vfmin
                  nm = nm + 1;
                  mats(nm) = m;
                end
              end
              mats = mats(1:nm);

              if nm == 1
                % Single material case
                obj.nmat_for_edge(k, j, i) = 1;
                dvol = sign_frac*vol_cell * abs(frac);
                dmass = dvol * mesh.rho_3dmat(k, j, i, mats(1));
                dener = dvol * mesh.ei_3dmat(k, j, i, mats(1));

                obj.matid_for_edge(loc_edge) = mats(1);
                obj.dvol_for_edge(loc_edge) = dvol;
                obj.dmass_for_edge(loc_edge) = dmass;
                obj.dener_for_edge(loc_edge) = dener;
                loc_edge = loc_edge + 1;
              else
                % Multiple material case: build polyhedron and advect
                xl_slab(dir) = xyz_edge;
                xr_slab(dir) = xyz_edge + abs(dist);
                xl_cell(dir) = xyz_edge;
                inward_norm(dir) = -1.0;
                slab_faceid = 2*dir; % 2, 4, 6
                
                ijk = [i, j, k];

                %% Build polyhedron
                [nnode_tot, coords_tot, nface_for_mpoly, ...
                  nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly] = ...
                  obj.build_mpoly3d(mesh, ijk);

                % Perform advection
                [nmat_adv, matid_adv, vol_adv, mass_adv, ener_adv] = ...
                  mat.advect3d(mesh, ijk, xl_cell, dx, nm, mats, ...
                  xl_slab, xr_slab, inward_norm, slab_faceid, ...
                  nnode_tot, coords_tot, nface_for_mpoly, ...
                  nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly);

                % Free temporary arrays
                coords_tot = [];
                nnode_for_face_ea_mpoly = [];
                nodelist_for_face_ea_mpoly = [];

                % Store results
                obj.nmat_for_edge(k, j, i) = nmat_adv;
                obj.matid_for_edge(loc_edge:loc_edge + nmat_adv - 1) = matid_adv(1:nmat_adv);
                obj.dvol_for_edge(loc_edge:loc_edge + nmat_adv - 1)  = sign_frac*vol_adv(1:nmat_adv);
                obj.dmass_for_edge(loc_edge:loc_edge + nmat_adv - 1) = sign_frac*mass_adv(1:nmat_adv);
                obj.dener_for_edge(loc_edge:loc_edge + nmat_adv - 1) = sign_frac*ener_adv(1:nmat_adv);
                loc_edge = loc_edge + nmat_adv;

              end
            end  % if abs(frac) > obj.small

            % if h5_checking
            %   global ncycle;
            %   aux = load_h5debug_3d('debug3d.h5',sprintf("dmass_for_edge(dir=%d)",(dir-1)),ncycle);
            %   c_dmass = squeeze(aux);
            %   if abs(obj.dmass_for_edge(loc_edge-1) - c_dmass(loc_edge-1))>1e-8
            %     fprintf(' Sign mismatch at loc_edge = %d\n', loc_edge-1);
            %     fprintf('vol_c = %.5e, vol_m = %.5e\n', c_dmass(loc_edge-1), obj.dmass_for_edge(loc_edge-1));
            %     keyboard  % ⛔ Pauses execution and opens the debugger
            %   end
            % end

          end % i
        end % j
      end % k
    end

    %%
    function update_materials_after_advection(obj, dir, nmat, vol_cell, sizes_edge, ncycle)
      % update_materials_after_advection - Update material quantities after advection.
      %
      % This function updates the material volume fractions, densities, and energies
      % in the computational domain after advection is performed in a given direction.
      %
      % Inputs:
      %   obj        - (c_advection) The advection object instance.
      %   dir        - (int) Direction of advection (1, 2, or 3 for x, y, z).
      %   nmat       - (int) Total number of materials in the simulation.
      %   vol_cell   - (double) Volume of a single cell (dx*dy*dz)
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

      global mesh mat h5_checking tol;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      ncell_ext = ncell + 2*nbdry;  % Add boundary cells
      ncell_bdry = ncell_ext - nbdry - 1;      % Define boundary limits
      loc_for_3dedge = zeros(sizes_edge(3), sizes_edge(2), sizes_edge(1));

      loc_edge = 1;
      for k = 1:sizes_edge(3)
        for j = 1:sizes_edge(2)
          for i = 1:sizes_edge(1)
            if obj.nmat_for_edge(k, j, i) > 0
              loc_for_3dedge(k, j, i) = loc_edge;
              loc_edge = loc_edge + obj.nmat_for_edge(k, j, i);
            else
              loc_for_3dedge(k, j, i) = -1;
            end
          end
        end
      end

      % Loop through all cells in the extended domain
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

          end %i
        end %j
      end %k

      %%
     
      if h5_checking
        aux = load_h5debug_4d('debug3d.h5',sprintf("vol_for_cell_after_advec(dir=%d)",(dir-1)),ncycle);
        assert(norm(obj.vol_for_cell-aux,'fro')<tol);
        printTextWithOK(sprintf('dir=%d - vol_for_cell_after_advec ',dir),55);

        aux = load_h5debug_4d('debug3d.h5',sprintf("mass_for_cell_after_advec(dir=%d)",(dir-1)),ncycle);
        assert(norm(obj.mass_for_cell-aux,'fro')<tol);
        printTextWithOK(sprintf('dir=%d - mass_for_cell_after_advec ',dir),55);

        aux = load_h5debug_4d('debug3d.h5',sprintf("ener_for_cell_after_advec(dir=%d)",(dir-1)),ncycle);
        assert(norm(obj.ener_for_cell-aux,'fro')<tol);
        printTextWithOK(sprintf('dir=%d - ener_for_cell_after_advec ',dir),55);

      end

      %% Remove small volume fractions
      for k = nbdry+1:ncell_bdry(3)+1
        for j = nbdry+1:ncell_bdry(2)+1
          for i = nbdry+1:ncell_bdry(1)+1
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

          end %i
        end
      end
      clear loc_for_3dedge;

      if h5_checking
        aux = load_h5debug_4d('debug3d.h5',sprintf("vol_for_cell_remove_frac(dir=%d)",(dir-1)),ncycle);
        assert(norm(obj.vol_for_cell-aux,'fro')<tol);
        % printTextWithOK(sprintf('dir=%d - vol_for_cell_remove_frac ',dir),55);

        aux = load_h5debug_4d('debug3d.h5',sprintf("mass_for_cell_remove_frac(dir=%d)",(dir-1)),ncycle);
        assert(norm(obj.mass_for_cell-aux,'fro')<tol);
        % printTextWithOK(sprintf('dir=%d - mass_for_cell_remove_frac ',dir),55);

        aux = load_h5debug_4d('debug3d.h5',sprintf("ener_for_cell_remove_frac(dir=%d)",(dir-1)),ncycle);
        assert(norm(obj.ener_for_cell-aux,'fro')<tol);
        % printTextWithOK(sprintf('dir=%d - ener_for_cell_after_remove_frac ',dir),55);

        printTextWithOK(sprintf('dir=%d - check cell variables after remove frac ',dir),55);
      end

    end

    %%
    function [nnode_tot, coords_tot, nnode_for_mpoly, nodes_for_mpoly] = ...
        build_mpoly2d(~, nbdry, xl_prob, dx, nmat, vf_2dmat, ijk, nm_this_cell)
      % build_mpoly2d Constructs material polygons for a 2D cell.
      %
      % Inputs:
      %   - nbdry: Number of boundary layers.
      %   - xl_prob: Lower bounds of the problem domain.
      %   - dx: Cell dimensions.
      %   - nmat: Number of materials.
      %   - vf_2dmat: Volume fractions for each material in each cell.
      %   - ijk: Cell indices [i, j].
      %   - nm_this_cell: Number of materials in the target cell.
      %
      % Outputs:
      %   - nnode_tot: Total number of nodes in the polygons.
      %   - coords_tot: Coordinates of the nodes.
      %   - nnode_for_mpoly: Number of nodes per material polygon.
      %   - nodes_for_mpoly: Indices of the nodes for each material polygon.

      global mat vof2d;
      geop = 1;
      dim = 2;

      % Calculate the lower corner coordinates of the current cell
      xl = xl_prob + (ijk - (nbdry+1)) .* dx(1:dim);
      icell = ijk(1);
      jcell = ijk(2);
      % vfs = vf_2dmat{jcell, icell};

      % Initialize nmat_2dsmesh, matid_2dsmesh, and vf_2dsmesh
      nmat_2dsmesh = zeros(3, 3);
      matid_2dsmesh = -ones(3, 3, nmat);
      vf_2dsmesh = zeros(3, 3, nmat);

      % Populate 2D stencil with material data
      for j0 = 1:3
        j = jcell + j0 - 2;  % Adjusted index for boundary
        for i0 = 1:3
          i = icell + i0 - 2;  % Adjusted index for boundary
          myvfs = vf_2dmat(j, i,:);
          nm = 0;
          for m = 1:nmat
            if myvfs(m) >= mat.vfmin
              nm = nm + 1;
              matid_2dsmesh(j0, i0, nm) = m;
              vf_2dsmesh(j0, i0, nm) = myvfs(m);
            end
          end
          nmat_2dsmesh(j0, i0) = nm;
        end
      end

      % Ensure the number of materials matches the target cell
      assert(nm_this_cell == nmat_2dsmesh(2, 2));

      % Initialize total number of nodes and material interfaces
      nnode_tot = 0;

      % Call reconstruct2d_nmat_pagosa for polygon reconstruction
      % [coords_tot, nnode_for_mpoly, nodes_for_mpoly] = vof2d.reconstruct2d_nmat_pagosa(geop, xl, dx, ...
      %   nmat_2dsmesh, matid_2dsmesh, vf_2dsmesh);

      [coords_tot, nnode_tot, nodes_for_minterface, ...
        nodes_for_mpoly, nnode_for_minterface, nnode_for_mpoly] ...
        = vof2d.reconstruct2d_nmat_pagosa(geop, xl, dx, ...
        nmat_2dsmesh, matid_2dsmesh, vf_2dsmesh);

    end

  end
end
