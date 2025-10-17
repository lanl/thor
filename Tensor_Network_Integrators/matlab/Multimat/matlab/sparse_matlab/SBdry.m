classdef SBdry < handle
  %   c_bdry  Boundary conditions

  properties
    nmat_for_cell = [];
    matid_for_cell = {};
    vol_for_cell = {};
    mass_for_cell = {};
    ener_for_cell = {};
  end

  methods

    %%%%%%%%%%%%%%% implemented & tested
    %%
    function bdry_cell_2d(~,mesh, btype_lower, btype_upper)
      % bdry_cell_2d Applies boundary conditions to 2D cells for transmitted boundaries.

      % Extract relevant properties
      ncell = mesh.ncell_prob;     % Number of cells in each direction
      nbdry = mesh.nbdry_prob;     % Number of boundary cells
      ncell_ext = ncell + 2 * nbdry; % Extended cell dimensions including boundaries
      nmat = mesh.nmat_mesh;       % Number of materials

      % xl edge without corners
      if btype_lower(1) == bdry_type.bdry_transmitted
        i0 = nbdry + 1;
        for j = nbdry+1:ncell(2)+nbdry
          for i = 1:nbdry
            % Initialize cell properties
            mesh.rho_for_2dcell(j, i) = 0.0;
            mesh.ei_for_2dcell(j, i) = 0.0;
            mesh.pres_for_2dcell(j, i) = 0.0;

            % Loop over materials
            for m = 1:nmat
              % Update material properties
              mesh.vf_2dmat(j, i, m) = mesh.vf_2dmat(j, i0, m);
              mesh.rho_2dmat(j, i, m) = mesh.rho_2dmat(j, i0, m);
              mesh.ei_2dmat(j, i, m) = mesh.ei_2dmat(j, i0, m);
              mesh.pres_2dmat(j, i, m) = mesh.pres_2dmat(j, i0, m);

              % Aggregate properties for the cell
              mesh.rho_for_2dcell(j, i) = mesh.rho_for_2dcell(j, i) + ...
                mesh.vf_2dmat(j, i, m) * mesh.rho_2dmat(j, i, m);
              mesh.ei_for_2dcell(j, i) = mesh.ei_for_2dcell(j, i) + ...
                mesh.vf_2dmat(j, i, m) * mesh.ei_2dmat(j, i, m);
              mesh.pres_for_2dcell(j, i) = mesh.pres_for_2dcell(j, i) + ...
                mesh.vf_2dmat(j, i, m) * mesh.pres_2dmat(j, i, m);
            end
          end
        end
      end

      % xr edge without corners
      if btype_upper(1) == bdry_type.bdry_transmitted
        i0 = ncell(1) + nbdry;
        for j = nbdry+1:ncell(2)+nbdry
          for i = 1:nbdry
            i1 = i0 + i;

            % Initialize cell properties
            mesh.rho_for_2dcell(j, i1) = 0.0;
            mesh.ei_for_2dcell(j, i1) = 0.0;
            mesh.pres_for_2dcell(j, i1) = 0.0;

            % Loop over materials
            for m = 1:nmat
              % Update material properties
              mesh.vf_2dmat(j, i1, m) = mesh.vf_2dmat(j, i0, m);
              mesh.rho_2dmat(j, i1, m) = mesh.rho_2dmat(j, i0, m);
              mesh.ei_2dmat(j, i1, m) = mesh.ei_2dmat(j, i0, m);
              mesh.pres_2dmat(j, i1, m) = mesh.pres_2dmat(j, i0, m);

              % Aggregate properties for the cell
              mesh.rho_for_2dcell(j, i1) = mesh.rho_for_2dcell(j, i1) + ...
                mesh.vf_2dmat(j, i1, m) * mesh.rho_2dmat(j, i1, m);
              mesh.ei_for_2dcell(j, i1) = mesh.ei_for_2dcell(j, i1) + ...
                mesh.vf_2dmat(j, i1, m) * mesh.ei_2dmat(j, i1, m);
              mesh.pres_for_2dcell(j, i1) = mesh.pres_for_2dcell(j, i1) + ...
                mesh.vf_2dmat(j, i1, m) * mesh.pres_2dmat(j, i1, m);
            end
          end
        end
      end

      % yl edge including corners
      if btype_lower(2) == bdry_type.bdry_transmitted
        j0 = nbdry + 1;
        for j = 1:nbdry
          for i = 1:ncell_ext(1)
            % Initialize cell properties
            mesh.rho_for_2dcell(j, i) = 0.0;
            mesh.ei_for_2dcell(j, i) = 0.0;
            mesh.pres_for_2dcell(j, i) = 0.0;

            % Loop over materials
            for m = 1:nmat
              % Update material properties
              mesh.vf_2dmat(j, i, m) = mesh.vf_2dmat(j0, i, m);
              mesh.rho_2dmat(j, i, m) = mesh.rho_2dmat(j0, i, m);
              mesh.ei_2dmat(j, i, m) = mesh.ei_2dmat(j0, i, m);
              mesh.pres_2dmat(j, i, m) = mesh.pres_2dmat(j0, i, m);

              % Aggregate properties for the cell
              mesh.rho_for_2dcell(j, i) = mesh.rho_for_2dcell(j, i) + ...
                mesh.vf_2dmat(j, i, m) * mesh.rho_2dmat(j, i, m);
              mesh.ei_for_2dcell(j, i) = mesh.ei_for_2dcell(j, i) + ...
                mesh.vf_2dmat(j, i, m) * mesh.ei_2dmat(j, i, m);
              mesh.pres_for_2dcell(j, i) = mesh.pres_for_2dcell(j, i) + ...
                mesh.vf_2dmat(j, i, m) * mesh.pres_2dmat(j, i, m);
            end
          end
        end
      end

      % yr edge including corners
      if btype_upper(2) == bdry_type.bdry_transmitted
        j0 = ncell(2) + nbdry;
        for j = 1:nbdry
          j1 = j0 + j;
          for i = 1:ncell_ext(1)
            % Initialize cell properties
            mesh.rho_for_2dcell(j1, i) = 0.0;
            mesh.ei_for_2dcell(j1, i) = 0.0;
            mesh.pres_for_2dcell(j1, i) = 0.0;

            % Loop over materials
            for m = 1:nmat
              % Update material properties
              mesh.vf_2dmat(j1, i, m) = mesh.vf_2dmat(j0, i, m);
              mesh.rho_2dmat(j1, i, m) = mesh.rho_2dmat(j0, i, m);
              mesh.ei_2dmat(j1, i, m) = mesh.ei_2dmat(j0, i, m);
              mesh.pres_2dmat(j1, i, m) = mesh.pres_2dmat(j0, i, m);

              % Aggregate properties for the cell
              mesh.rho_for_2dcell(j1, i) = mesh.rho_for_2dcell(j1, i) + ...
                mesh.vf_2dmat(j1, i, m) * mesh.rho_2dmat(j1, i, m);
              mesh.ei_for_2dcell(j1, i) = mesh.ei_for_2dcell(j1, i) + ...
                mesh.vf_2dmat(j1, i, m) * mesh.ei_2dmat(j1, i, m);
              mesh.pres_for_2dcell(j1, i) = mesh.pres_for_2dcell(j1, i) + ...
                mesh.vf_2dmat(j1, i, m) * mesh.pres_2dmat(j1, i, m);
            end
          end
        end
      end
    end

    function bdry_node_2d(~, mesh, btype_lower, btype_upper)
      % bdry_node_2d Applies boundary conditions to 2D node velocities.
      %
      % This method handles transmitted boundary conditions for node velocities
      % on the edges of a 2D grid, updating the velocity arrays for boundary
      % nodes without modifying the interior.
      %
      % Inputs:
      %   mesh        - (c_mesh) The mesh object containing grid and velocity data.
      %   btype_lower - (array of int) Boundary type for the lower boundaries (xl, yl).
      %   btype_upper - (array of int) Boundary type for the upper boundaries (xr, yr).
      %
      % Uses:
      %   mesh.ncell_prob      - (array of int) Number of cells in each spatial dimension.
      %   mesh.nbdry_prob      - (int) Number of boundary cells.
      %   mesh.vel_for_2dnode  - (3D array) Velocities at each 2D node.
      %   mesh.vav_for_2dnode - (3D array) Averaged velocities at each 2D node.
      %
      % Modifies:
      %   mesh.vel_for_2dnode  - Updates the velocities at the boundary nodes.
      %   mesh.vav_for_2dnode - Updates the averaged velocities at the boundary nodes.
      %
      % Original C declaration:
      % void bdry_node_2d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                   double ***vel_for_2dnode, double ***vav_for_2dnode);

      dim = 2;
      ncell = mesh.ncell_prob;   % Number of cells in each direction
      nbdry = mesh.nbdry_prob;   % Number of boundary cells
      nnode = ncell + 1;
      nnode_ext = ncell + 2 * nbdry + 1;  % Extended grid size (including boundaries)

      % xl edge without corners
      if btype_lower(1) == bdry_type.bdry_transmitted
        n0 = nbdry + 1;
        for j = nbdry + 1:nnode(2) + nbdry
          for i = 1:nbdry
            for dir = 1:dim
              mesh.vel_for_2dnode(j, i, dir) = mesh.vel_for_2dnode(j, n0, dir);
              mesh.vav_for_2dnode(j, i, dir) = mesh.vav_for_2dnode(j, n0, dir);
            end
          end
        end
      end

      % xr edge without corners
      if btype_upper(1) == bdry_type.bdry_transmitted
        n0 = nnode(1) + nbdry;
        for j = nbdry + 1:nnode(2) + nbdry
          for i = 1:nbdry
            n1 = n0 + i;
            for dir = 1:dim
              mesh.vel_for_2dnode(j, n1, dir) = mesh.vel_for_2dnode(j, n0, dir);
              mesh.vav_for_2dnode(j, n1, dir) = mesh.vav_for_2dnode(j, n0, dir);
            end
          end
        end
      end

      % yl edge including corners
      if btype_lower(2) == bdry_type.bdry_transmitted
        n0 = nbdry + 1;
        for j = 1:nbdry
          for i = 1:nnode_ext(1)
            for dir = 1:dim
              mesh.vel_for_2dnode(j, i, dir) = mesh.vel_for_2dnode(n0, i, dir);
              mesh.vav_for_2dnode(j, i, dir) = mesh.vav_for_2dnode(n0, i, dir);
            end
          end
        end
      end

      % yr edge including corners
      if btype_upper(2) == bdry_type.bdry_transmitted
        n0 = nnode(2) + nbdry;
        for j = 1:nbdry
          n1 = n0 + j;
          for i = 1:nnode_ext(1)
            for dir = 1:dim
              mesh.vel_for_2dnode(n1, i, dir) = mesh.vel_for_2dnode(n0, i, dir);
              mesh.vav_for_2dnode(n1, i, dir) = mesh.vav_for_2dnode(n0, i, dir);
            end
          end
        end
      end
    end

    function bdry_cell_1var_2d(obj, mesh, btype_lower, btype_upper)
      % bdry_cell_1var_2d Applies boundary conditions for a 2D cell array of materials.
      %
      % This function applies transmitted boundary conditions for the number of
      % materials in a 2D mesh, handling the lower and upper boundaries for both
      % the x and y directions, with specific rules for edges and corners.
      %
      % Inputs:
      %   mesh          - (c_mesh) The mesh object containing the number of cells and boundary data.
      %   btype_lower   - (array of int) Boundary type for the lower boundaries (xl, yl).
      %   btype_upper   - (array of int) Boundary type for the upper boundaries (xr, yr).
      %
      % Outputs:
      %   None. The function modifies the number of materials in each cell at the boundaries.
      %
      % Uses:
      %   mesh.ncell_prob      - (array of int) Number of cells in each direction [size: 1x2].
      %   mesh.nbdry_prob      - (int) Number of boundary layers.
      %
      % Modifies:
      %   bdry.nmat_for_cell - (2D array of int) Updates the number of materials in cells
      %                          at the boundaries (xl, xr, yl, yr).
      %
      % Original C declaration:
      % void bdry_cell_1var_2d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                        int **nmat_for_2dcell)

      dim = 2;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      ncell_ext = ncell + 2 * nbdry;  % Extended number of cells including boundary layers

      % xl edge without corners
      if btype_lower(1) == bdry_type.bdry_transmitted
        i0 = nbdry + 1;
        for j = nbdry+1:ncell(2) + nbdry
          for i = 1:nbdry
            obj.nmat_for_cell(j, i) = obj.nmat_for_cell(j, i0);
          end
        end
      end

      % xr edge without corners
      if btype_upper(1) == bdry_type.bdry_transmitted
        i0 = ncell(1) + nbdry;
        for j = nbdry+1:ncell(2) + nbdry
          for i = 1:nbdry
            i1 = i0 + i;
            obj.nmat_for_cell(j, i1) = obj.nmat_for_cell(j, i0);
          end
        end
      end

      % yl edge including corners
      if btype_lower(2) == bdry_type.bdry_transmitted
        j0 = nbdry + 1;
        for j = 1:nbdry
          for i = 1:ncell_ext(1)
            obj.nmat_for_cell(j, i) = obj.nmat_for_cell(j0, i);
          end
        end
      end

      % yr edge including corners
      if btype_upper(2) == bdry_type.bdry_transmitted
        j0 = ncell(2) + nbdry;
        for j = 1:nbdry
          j1 = j0 + j;
          for i = 1:ncell_ext(1)
            obj.nmat_for_cell(j1, i) = obj.nmat_for_cell(j0, i);
          end
        end
      end
    end

    function bdry_cell_ragged_2d(obj, mesh, btype_lower, btype_upper)
      % bdry_cell_ragged_2d Applies boundary conditions for a 2D ragged cell array of materials.
      %
      % This method handles transmitted boundary conditions for the number of
      % materials and their properties (volume, mass, energy) in a 2D mesh,
      % updating boundary cells by copying data from inner cells.
      %
      % Inputs:
      %   obj           - (c_bdry) The current boundary object instance.
      %   mesh          - (c_mesh) The mesh object containing the cell data.
      %   btype_lower   - (array of int) Boundary type for the lower boundaries (xl, yl).
      %   btype_upper   - (array of int) Boundary type for the upper boundaries (xr, yr).
      %
      % Outputs:
      %   None. The function modifies the material properties at the boundaries.
      %
      % Uses:
      %   mesh.ncell_prob       - (array of int) Number of cells in each direction [size: 1x2].
      %   mesh.nbdry_prob       - (int) Number of boundary layers.
      %
      % Modifies:
      %   obj.nmat_for_cell     - (2D array of double) Number of materials for each cell.
      %   obj.matid_for_cell    - (cell array of int arrays) Material IDs for each cell.
      %   obj.vol_for_cell      - (cell array of double arrays) Volume for each material in each cell.
      %   obj.mass_for_cell     - (cell array of double arrays) Mass for each material in each cell.
      %   obj.ener_for_cell     - (cell array of double arrays) Energy for each material in each cell.
      %
      % Original C declaration:
      % void bdry_cell_ragged_2d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                          int **nmat_for_cell, int ***matid_for_cell, double ***vol_for_cell,
      %                          double ***mass_for_cell, double ***ener_for_cell)

      dim = 2;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      ncell_ext = ncell + 2 * nbdry;  % Extended number of cells including boundary layers

      % xl edge without corners
      if btype_lower(1) == bdry_type.bdry_transmitted
        i0 = nbdry + 1;
        for j = nbdry+1:ncell(2) + nbdry
          for i = 1:nbdry
            obj.nmat_for_cell(j, i) = obj.nmat_for_cell(j, i0);
            nm = obj.nmat_for_cell(j, i0);
            obj.matid_for_cell{j, i}(1:nm) = obj.matid_for_cell{j, i0}(1:nm);
            obj.vol_for_cell{j, i}(1:nm) = obj.vol_for_cell{j, i0}(1:nm);
            obj.mass_for_cell{j, i}(1:nm) = obj.mass_for_cell{j, i0}(1:nm);
            obj.ener_for_cell{j, i}(1:nm) = obj.ener_for_cell{j, i0}(1:nm);
          end
        end
      end

      % xr edge without corners
      if btype_upper(1) == bdry_type.bdry_transmitted
        i0 = ncell(1) + nbdry;
        for j = nbdry+1:ncell(2) + nbdry
          for i = 1:nbdry
            i1 = i0 + i;
            obj.nmat_for_cell(j, i1) = obj.nmat_for_cell(j, i0);
            nm = obj.nmat_for_cell(j, i0);
            obj.matid_for_cell{j, i1}(1:nm) = obj.matid_for_cell{j, i0}(1:nm);
            obj.vol_for_cell{j, i1}(1:nm) = obj.vol_for_cell{j, i0}(1:nm);
            obj.mass_for_cell{j, i1}(1:nm) = obj.mass_for_cell{j, i0}(1:nm);
            obj.ener_for_cell{j, i1}(1:nm) = obj.ener_for_cell{j, i0}(1:nm);
          end
        end
      end

      % yl edge including corners
      if btype_lower(2) == bdry_type.bdry_transmitted
        j0 = nbdry + 1;
        for j = 1:nbdry
          for i = 1:ncell_ext(1)
            obj.nmat_for_cell(j, i) = obj.nmat_for_cell(j0, i);
            nm = obj.nmat_for_cell(j0, i);
            obj.matid_for_cell{j, i}(1:nm) = obj.matid_for_cell{j0, i}(1:nm);
            obj.vol_for_cell{j, i}(1:nm) = obj.vol_for_cell{j0, i}(1:nm);
            obj.mass_for_cell{j, i}(1:nm) = obj.mass_for_cell{j0, i}(1:nm);
            obj.ener_for_cell{j, i}(1:nm) = obj.ener_for_cell{j0, i}(1:nm);
          end
        end
      end

      % yr edge including corners
      if btype_upper(2) == bdry_type.bdry_transmitted
        j0 = ncell(2) + nbdry;
        for j = 1:nbdry
          j1 = j0 + j;
          for i = 1:ncell_ext(1)
            obj.nmat_for_cell(j1, i) = obj.nmat_for_cell(j0, i);
            nm = obj.nmat_for_cell(j0, i);
            obj.matid_for_cell{j1, i}(1:nm) = obj.matid_for_cell{j0, i}(1:nm);
            obj.vol_for_cell{j1, i}(1:nm) = obj.vol_for_cell{j0, i}(1:nm);
            obj.mass_for_cell{j1, i}(1:nm) = obj.mass_for_cell{j0, i}(1:nm);
            obj.ener_for_cell{j1, i}(1:nm) = obj.ener_for_cell{j0, i}(1:nm);
          end
        end
      end
    end

    %%%%%%%%%%%%%%% implemented but not tested

    function cell_bdry_condition(obj, mesh, btype_lower, btype_upper)
      % cell_bdry_condition: simple switch

      if mesh.dim_prob == 2
        obj.bdry_cell_2d(mesh, btype_lower, btype_upper);
        obj.bdry_node_2d(mesh, btype_lower, btype_upper);

      elseif mesh.dim_prob == 3
        obj.bdry_cell_3d(mesh, btype_lower, btype_upper);
        obj.bdry_node_3d(mesh, btype_lower, btype_upper);

      end
    end

    function mixcell_bdry_condition(~, mesh, mat, btype_lower, btype_upper)
      % mixcell_bdry_condition  Applies boundary conditions to mixed cells in a 2D or 3D mesh.
      %
      % This method adjusts the mixed cells at the boundaries of the mesh based on the
      % provided boundary types. It handles boundary-specific adjustments in both 2D
      % and 3D cases and updates material properties such as densities, pressures,
      % internal energies, and volume fractions for the mixed cells.
      %
      % Inputs:
      %   obj          - (c_bdry) The boundary object instance.
      %   mesh         - (c_mesh) The mesh object instance containing grid information.
      %   mat          - (c_mat) The material object instance that stores mixed cell properties.
      %   btype_lower  - (array of int) Boundary type for the lower boundary in each direction.
      %   btype_upper  - (array of int) Boundary type for the upper boundary in each direction.
      %
      % Outputs:
      %   None. The function modifies mixed cell properties in the `mesh` and `mat` objects.
      %
      % Uses:
      %   mesh.dim           - (int) Dimensionality of the problem (2D or 3D).
      %   mesh.ncell_prob    - (array of int) Number of cells in each spatial dimension.
      %   mesh.nbdry_prob    - (int) Number of boundary cells.
      %   mat.nmixcell_int   - (int) Number of mixed cells in interior cells.
      %
      % Modifies:
      %   mat.nmat_in_mixcell   - (array of int) Number of materials in each mixed cell.
      %   mat.ijk_in_mixcell    - (cell of int arrays) Indices of mixed cells in the mesh.
      %   mat.matids_in_mixcell - (cell of int arrays) Material IDs in each mixed cell.
      %   mat.vf_in_mixcell     - (cell of double arrays) Volume fractions for materials in each mixed cell.
      %   mat.rho_in_mixcell    - (cell of double arrays) Densities of materials in each mixed cell.
      %   mat.pres_in_mixcell   - (cell of double arrays) Pressures of materials in each mixed cell.
      %   mat.ei_in_mixcell     - (cell of double arrays) Internal energies of materials in each mixed cell.
      %
      % Original C declaration:
      %   void mixcell_bdry_condition(int dim, int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                               int *nmixcell, int nmixcell_int, int **pnmat_in_mixcell, int ***pijk_in_mixcell,
      %                               int ***pmatids_in_mixcell, double ***pvf_in_mixcell, double ***prho_in_mixcell,
      %                               double ***ppres_in_mixcell, double ***pei_in_mixcell);

      % Initialize variables from mesh properties
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;

      % Calculate boundary cell limits
      ncell_bdry = ncell + nbdry - 1;

      % Maximum number of mixed cells
      nmix_max = mat.nmixcell_int;

      if dim == 2
        nmix_max = nmix_max + (4 * nbdry^2 + 2 * nbdry * (ncell(1) + ncell(2)));
      elseif dim == 3
        nmix_max = nmix_max + (8 * nbdry^3 + ...
          4 * nbdry^2 * sum(ncell) + ...
          2 * nbdry * (ncell(1) * ncell(2) + ...
          ncell(2) * ncell(3) + ...
          ncell(1) * ncell(3)));
      end

      % Reallocate memory for boundary cells and mixed cells
      mat.nmat_in_mixcell = [mat.nmat_in_mixcell; zeros(nmix_max - mat.nmixcell_int, 1)];
      mat.ijk_in_mixcell = [mat.ijk_in_mixcell; zeros(nmix_max, dim)];
      mat.matids_in_mixcell = [mat.matids_in_mixcell; cell(nmix_max, 1)];
      mat.vf_in_mixcell = [mat.vf_in_mixcell; cell(nmix_max, 1)];
      mat.rho_in_mixcell = [mat.rho_in_mixcell; cell(nmix_max, 1)];
      mat.pres_in_mixcell = [mat.pres_in_mixcell; cell(nmix_max, 1)];
      mat.ei_in_mixcell = [mat.ei_in_mixcell; cell(nmix_max, 1)];

      nmixcell_bdry = 0;
      nm_bdry = 0;
      offset = mat.nmixcell_int;

      % Initialize a variable to keep track of the additional boundary cells
      list_offset = mat.nmixcell_int;

      if dim == 2
        for mx = 1:mat.nmixcell_int
          i = mat.ijk_in_mixcell(mx, 1);
          j = mat.ijk_in_mixcell(mx, 2);

          if ((i == nbdry || i == ncell_bdry(1)) && (j == nbdry || j == ncell_bdry(2)))
            % Handle corner cells
            nmix = (nbdry + 1)^2 - 1;
            for idx = 1:nmix
              mat.nmat_in_mixcell(list_offset + idx - 1) = mat.nmat_in_mixcell(mx);
            end
            nmixcell_bdry = nmixcell_bdry + nmix;
            nm_bdry = nm_bdry + nmix * mat.nmat_in_mixcell(mx);
            list_offset = list_offset + nmix;

          elseif (i == nbdry || i == ncell_bdry(1) || j == nbdry || j == ncell_bdry(2))
            % Handle edge cells
            nmix = nbdry;
            for idx = 1:nmix
              mat.nmat_in_mixcell(list_offset + idx - 1) = mat.nmat_in_mixcell(mx);
            end
            nmixcell_bdry = nmixcell_bdry + nmix;
            nm_bdry = nm_bdry + nmix * mat.nmat_in_mixcell(mx);
            list_offset = list_offset + nmix;
          end
        end
      elseif dim == 3
        ntime_c = 0;
        ntime_f = 0;

        for mx = 1:mat.nmixcell_int
          i = mat.ijk_in_mixcell(mx, 1);
          j = mat.ijk_in_mixcell(mx, 2);
          k = mat.ijk_in_mixcell(mx, 3);

          % Handle corner cells
          if ((i == nbdry || i == ncell_bdry(1)) && ...
              (j == nbdry || j == ncell_bdry(2)) && ...
              (k == nbdry || k == ncell_bdry(3)))
            nmix = (nbdry + 1)^3 - 1;
            for idx = 1:nmix
              mat.nmat_in_mixcell(list_offset + idx - 1) = mat.nmat_in_mixcell(mx);
            end
            nmixcell_bdry = nmixcell_bdry + nmix;
            nm_bdry = nm_bdry + nmix * mat.nmat_in_mixcell(mx);
            list_offset = list_offset + nmix;

            % Handle z-edge cells
          elseif ((i == nbdry || i == ncell_bdry(1)) && ...
              (j == nbdry || j == ncell_bdry(2)))
            nmix = (nbdry + 1)^2 - 1;
            for idx = 1:nmix
              mat.nmat_in_mixcell(list_offset + idx - 1) = mat.nmat_in_mixcell(mx);
            end
            nmixcell_bdry = nmixcell_bdry + nmix;
            nm_bdry = nm_bdry + nmix * mat.nmat_in_mixcell(mx);
            list_offset = list_offset + nmix;

            % Handle y-edge cells
          elseif ((i == nbdry || i == ncell_bdry(1)) && ...
              (k == nbdry || k == ncell_bdry(3)))
            nmix = (nbdry + 1)^2 - 1;
            for idx = 1:nmix
              mat.nmat_in_mixcell(list_offset + idx - 1) = mat.nmat_in_mixcell(mx);
            end
            nmixcell_bdry = nmixcell_bdry + nmix;
            nm_bdry = nm_bdry + nmix * mat.nmat_in_mixcell(mx);
            list_offset = list_offset + nmix;

            % Handle x-edge cells
          elseif ((j == nbdry || j == ncell_bdry(2)) && ...
              (k == nbdry || k == ncell_bdry(3)))
            ntime_c = ntime_c + 1;
            nmix = (nbdry + 1)^2 - 1;
            for idx = 1:nmix
              mat.nmat_in_mixcell(list_offset + idx - 1) = mat.nmat_in_mixcell(mx);
            end
            nmixcell_bdry = nmixcell_bdry + nmix;
            nm_bdry = nm_bdry + nmix * mat.nmat_in_mixcell(mx);
            list_offset = list_offset + nmix;

            % Handle face cells (xl, xr, yl, yr, zl, zr)
          elseif ((i == nbdry || i == ncell_bdry(1)) || ...
              (j == nbdry || j == ncell_bdry(2)) || ...
              (k == nbdry || k == ncell_bdry(3)))
            ntime_f = ntime_f + 1;
            nmixcell_bdry = nmixcell_bdry + nbdry;
            nm_bdry = nm_bdry + nbdry * mat.nmat_in_mixcell(mx);
            list_offset = list_offset + nbdry;
          end
        end
      end

      if nmixcell_bdry == 0
        % If no boundary mixed cells, set nmixcell and return
        mat.nmixcell = mat.nmixcell_int;
        return;
      else
        % Sum the number of materials in mixed cells
        nm_sum = sum(mat.nmat_in_mixcell(1:mat.nmixcell_int));

        % Update nmixcell with boundary mixed cells
        mat.nmixcell = mat.nmixcell_int + nmixcell_bdry;
      end

      % Reallocate arrays to account for boundary mixed cells
      mat.nmat_in_mixcell = [mat.nmat_in_mixcell; zeros(mat.nmixcell - mat.nmixcell_int, 1)];
      mat.ijk_in_mixcell = [mat.ijk_in_mixcell; zeros(mat.nmixcell - mat.nmixcell_int, dim)];
      mat.matids_in_mixcell = [mat.matids_in_mixcell; cell(mat.nmixcell - mat.nmixcell_int, 1)];
      mat.vf_in_mixcell = [mat.vf_in_mixcell; cell(mat.nmixcell - mat.nmixcell_int, 1)];
      mat.rho_in_mixcell = [mat.rho_in_mixcell; cell(mat.nmixcell - mat.nmixcell_int, 1)];
      mat.pres_in_mixcell = [mat.pres_in_mixcell; cell(mat.nmixcell - mat.nmixcell_int, 1)];
      mat.ei_in_mixcell = [mat.ei_in_mixcell; cell(mat.nmixcell - mat.nmixcell_int, 1)];

      % Populate new allocations for boundary cells, from nmixcell_int to nmixcell
      for mx = mat.nmixcell_int + 1:mat.nmixcell
        nm = mat.nmat_in_mixcell(mx);  % Number of materials for this mixed cell

        % Allocate 1D arrays for each material property, based on nm
        mat.matids_in_mixcell{mx} = zeros(1, nm);
        mat.vf_in_mixcell{mx} = zeros(1, nm);
        mat.rho_in_mixcell{mx} = zeros(1, nm);
        mat.pres_in_mixcell{mx} = zeros(1, nm);
        mat.ei_in_mixcell{mx} = zeros(1, nm);
      end

      % Initialize ghost counters
      nghost_ylzl = 0;
      nghost_yrzl = 0;
      nghost_ylzr = 0;
      nghost_yrzr = 0;

      nghost_yl = 0;
      nghost_yr = 0;
      nghost_zl = 0;
      nghost_zr = 0;

      % Initialize the offset
      offset = mat.nmixcell_int;

      if dim == 2
        for mx = 1:mat.nmixcell_int
          i = mat.ijk_in_mixcell(mx, 1);
          j = mat.ijk_in_mixcell(mx, 2);

          if ((i == nbdry || i == ncell_bdry(1)) && (j == nbdry || j == ncell_bdry(2)))
            % Handle corner cells
            nmix = (nbdry + 1)^2 - 1;  % Exclude the corner cell itself
            nm = mat.nmat_in_mixcell(mx);
            mymx = offset;

            % xl yl corner
            if i == nbdry && j == nbdry
              for jb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && jb == 0, continue; end  % Exclude the corner cell itself
                  mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, nbdry - jb];
                  mymx = mymx + 1;
                end
              end
              % xl yr corner
            elseif i == nbdry && j == ncell_bdry(2)
              for jb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && jb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, ncell_bdry(2) + jb];
                  mymx = mymx + 1;
                end
              end
              % xr yl corner
            elseif i == ncell_bdry(1) && j == nbdry
              for jb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && jb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, nbdry - jb];
                  mymx = mymx + 1;
                end
              end
              % xr yr corner
            elseif i == ncell_bdry(1) && j == ncell_bdry(2)
              for jb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && jb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, ncell_bdry(2) + jb];
                  mymx = mymx + 1;
                end
              end
            end

            % Copy material properties for mixed cells
            for idx = 1:nmix
              mymx = offset + idx;
              for m = 1:nm
                mat.matids_in_mixcell{mymx}(m) = mat.matids_in_mixcell{mx}(m);
                mat.vf_in_mixcell{mymx}(m) = mat.vf_in_mixcell{mx}(m);
                mat.rho_in_mixcell{mymx}(m) = mat.rho_in_mixcell{mx}(m);
                mat.pres_in_mixcell{mymx}(m) = mat.pres_in_mixcell{mx}(m);
                mat.ei_in_mixcell{mymx}(m) = mat.ei_in_mixcell{mx}(m);
              end
            end

            % Update the offset for the next corner
            offset = offset + nmix;

          elseif (i == nbdry || i == ncell_bdry(1) || j == nbdry || j == ncell_bdry(2))  % edge
            nm = mat.nmat_in_mixcell(mx);
            mymx = offset;

            % Handle xl (left edge)
            if i == nbdry
              for ib = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [nbdry - 1 - ib, j];
                mymx = mymx + 1;
              end
              % Handle xr (right edge)
            elseif i == ncell_bdry(1)
              for ib = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + 1 + ib, j];
                mymx = mymx + 1;
              end
              % Handle yl (bottom edge)
            elseif j == nbdry
              for jb = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [i, nbdry - 1 - jb];
                mymx = mymx + 1;
              end
              % Handle yr (top edge)
            elseif j == ncell_bdry(2)
              for jb = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [i, ncell_bdry(2) + 1 + jb];
                mymx = mymx + 1;
              end
            end

            % Copy material properties for edge cells
            for idx = 1:nbdry
              mymx = offset + idx - 1;
              for m = 1:nm
                mat.matids_in_mixcell{mymx}(m) = mat.matids_in_mixcell{mx}(m);
                mat.vf_in_mixcell{mymx}(m) = mat.vf_in_mixcell{mx}(m);
                mat.rho_in_mixcell{mymx}(m) = mat.rho_in_mixcell{mx}(m);
                mat.pres_in_mixcell{mymx}(m) = mat.pres_in_mixcell{mx}(m);
                mat.ei_in_mixcell{mymx}(m) = mat.ei_in_mixcell{mx}(m);
              end
            end

            % Update the offset for the next edge
            offset = offset + nbdry;
          end

        end % for loop

      elseif dim == 3
        for mx = 1:mat.nmixcell_int
          i = mat.ijk_in_mixcell(mx, 1);
          j = mat.ijk_in_mixcell(mx, 2);
          k = mat.ijk_in_mixcell(mx, 3);

          if ((i == nbdry || i == ncell_bdry(1)) && ...
              (j == nbdry || j == ncell_bdry(2)) && ...
              (k == nbdry || k == ncell_bdry(3)))
            % Handle corner cells
            nmix = (nbdry + 1)^3 - 1;  % Exclude the corner cell itself
            nm = mat.nmat_in_mixcell(mx);
            mymx = offset;

            % Handle the xl yl zl corner
            if i == nbdry && j == nbdry && k == nbdry
              for kb = 0:nbdry
                for jb = 0:nbdry
                  for ib = 0:nbdry
                    if ib == 0 && jb == 0 && kb == 0, continue; end
                    mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, nbdry - jb, nbdry - kb];
                    mymx = mymx + 1;
                  end
                end
              end

              % Handle the xr yl zl corner
            elseif i == ncell_bdry(1) && j == nbdry && k == nbdry
              for kb = 0:nbdry
                for jb = 0:nbdry
                  for ib = 0:nbdry
                    if ib == 0 && jb == 0 && kb == 0, continue; end
                    mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, nbdry - jb, nbdry - kb];
                    mymx = mymx + 1;
                  end
                end
              end

              % Handle the xr yr zl corner
            elseif i == ncell_bdry(1) && j == ncell_bdry(2) && k == nbdry
              for kb = 0:nbdry
                for jb = 0:nbdry
                  for ib = 0:nbdry
                    if ib == 0 && jb == 0 && kb == 0, continue; end
                    mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, ncell_bdry(2) + jb, nbdry - kb];
                    mymx = mymx + 1;
                  end
                end
              end

              % Handle the xl yr zl corner
            elseif i == nbdry && j == ncell_bdry(2) && k == nbdry
              for kb = 0:nbdry
                for jb = 0:nbdry
                  for ib = 0:nbdry
                    if ib == 0 && jb == 0 && kb == 0, continue; end
                    mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, ncell_bdry(2) + jb, nbdry - kb];
                    mymx = mymx + 1;
                  end
                end
              end

              % Handle the xl yl zr corner
            elseif i == nbdry && j == nbdry && k == ncell_bdry(3)
              for kb = 0:nbdry
                for jb = 0:nbdry
                  for ib = 0:nbdry
                    if ib == 0 && jb == 0 && kb == 0, continue; end
                    mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, nbdry - jb, ncell_bdry(3) + kb];
                    mymx = mymx + 1;
                  end
                end
              end

              % Handle the xr yl zr corner
            elseif i == ncell_bdry(1) && j == nbdry && k == ncell_bdry(3)
              for kb = 0:nbdry
                for jb = 0:nbdry
                  for ib = 0:nbdry
                    if ib == 0 && jb == 0 && kb == 0, continue; end
                    mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, nbdry - jb, ncell_bdry(3) + kb];
                    mymx = mymx + 1;
                  end
                end
              end

              % Handle the xr yr zr corner
            elseif i == ncell_bdry(1) && j == ncell_bdry(2) && k == ncell_bdry(3)
              for kb = 0:nbdry
                for jb = 0:nbdry
                  for ib = 0:nbdry
                    if ib == 0 && jb == 0 && kb == 0, continue; end
                    mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, ncell_bdry(2) + jb, ncell_bdry(3) + kb];
                    mymx = mymx + 1;
                  end
                end
              end

              % Handle the xl yr zr corner
            elseif i == nbdry && j == ncell_bdry(2) && k == ncell_bdry(3)
              for kb = 0:nbdry
                for jb = 0:nbdry
                  for ib = 0:nbdry
                    if ib == 0 && jb == 0 && kb == 0, continue; end
                    mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, ncell_bdry(2) + jb, ncell_bdry(3) + kb];
                    mymx = mymx + 1;
                  end
                end
              end
            end

            % Copy material properties for corner cells
            for idx = 1:nmix
              mymx = offset + idx;
              for m = 1:nm
                mat.matids_in_mixcell{mymx}(m) = mat.matids_in_mixcell{mx}(m);
                mat.vf_in_mixcell{mymx}(m) = mat.vf_in_mixcell{mx}(m);
                mat.rho_in_mixcell{mymx}(m) = mat.rho_in_mixcell{mx}(m);
                mat.pres_in_mixcell{mymx}(m) = mat.pres_in_mixcell{mx}(m);
                mat.ei_in_mixcell{mymx}(m) = mat.ei_in_mixcell{mx}(m);
              end
            end

            % Update the offset for the next corner
            offset = offset + nmix;

          elseif (i == nbdry || i == ncell_bdry(1)) && (j == nbdry || j == ncell_bdry(2))  % z-edge
            nmix = (nbdry + 1)^2 - 1;  % Exclude the edge cell itself
            nm = mat.nmat_in_mixcell(mx);
            mymx = offset;

            % Handle xlyl edge
            if i == nbdry && j == nbdry
              for jb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && jb == 0, continue; end  % Exclude the edge cell itself
                  mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, nbdry - jb, k];
                  mymx = mymx + 1;
                end
              end

              % Handle xryl edge
            elseif i == ncell_bdry(1) && j == nbdry
              for jb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && jb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, nbdry - jb, k];
                  mymx = mymx + 1;
                end
              end

              % Handle xlyr edge
            elseif i == nbdry && j == ncell_bdry(2)
              for jb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && jb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, ncell_bdry(2) + jb, k];
                  mymx = mymx + 1;
                end
              end

              % Handle xryr edge
            elseif i == ncell_bdry(1) && j == ncell_bdry(2)
              for jb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && jb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, ncell_bdry(2) + jb, k];
                  mymx = mymx + 1;
                end
              end
            end

            % Copy material properties for edge cells
            for idx = 1:nmix
              mymx = offset + idx;
              for m = 1:nm
                mat.matids_in_mixcell{mymx}(m) = mat.matids_in_mixcell{mx}(m);
                mat.vf_in_mixcell{mymx}(m) = mat.vf_in_mixcell{mx}(m);
                mat.rho_in_mixcell{mymx}(m) = mat.rho_in_mixcell{mx}(m);
                mat.pres_in_mixcell{mymx}(m) = mat.pres_in_mixcell{mx}(m);
                mat.ei_in_mixcell{mymx}(m) = mat.ei_in_mixcell{mx}(m);
              end
            end

            % Update the offset for the next edge
            offset = offset + nmix;

          elseif (i == nbdry || i == ncell_bdry(1)) && (k == nbdry || k == ncell_bdry(3))  % y-edge
            nmix = (nbdry + 1)^2 - 1;  % Exclude the edge cell itself
            nm = mat.nmat_in_mixcell(mx);
            mymx = offset;

            % Handle xlzl edge
            if i == nbdry && k == nbdry
              for kb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && kb == 0, continue; end  % Exclude the edge cell itself
                  mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, j, nbdry - kb];
                  mymx = mymx + 1;
                end
              end

              % Handle xrzl edge
            elseif i == ncell_bdry(1) && k == nbdry
              for kb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && kb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, j, nbdry - kb];
                  mymx = mymx + 1;
                end
              end

              % Handle xlzr edge
            elseif i == nbdry && k == ncell_bdry(3)
              for kb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && kb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [nbdry - ib, j, ncell_bdry(3) + kb];
                  mymx = mymx + 1;
                end
              end

              % Handle xrzr edge
            elseif i == ncell_bdry(1) && k == ncell_bdry(3)
              for kb = 0:nbdry
                for ib = 0:nbdry
                  if ib == 0 && kb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + ib, j, ncell_bdry(3) + kb];
                  mymx = mymx + 1;
                end
              end
            end

            % Copy material properties for edge cells
            for idx = 1:nmix
              mymx = offset + idx;
              for m = 1:nm
                mat.matids_in_mixcell{mymx}(m) = mat.matids_in_mixcell{mx}(m);
                mat.vf_in_mixcell{mymx}(m) = mat.vf_in_mixcell{mx}(m);
                mat.rho_in_mixcell{mymx}(m) = mat.rho_in_mixcell{mx}(m);
                mat.pres_in_mixcell{mymx}(m) = mat.pres_in_mixcell{mx}(m);
                mat.ei_in_mixcell{mymx}(m) = mat.ei_in_mixcell{mx}(m);
              end
            end

            % Update the offset for the next edge
            offset = offset + nmix;

          elseif (j == nbdry || j == ncell_bdry(2)) && (k == nbdry || k == ncell_bdry(3))  % x-edge
            nmix = (nbdry + 1)^2 - 1;  % Exclude the edge cell itself
            nm = mat.nmat_in_mixcell(mx);
            mymx = offset;

            % Handle ylzl edge
            if j == nbdry && k == nbdry
              for kb = 0:nbdry
                for jb = 0:nbdry
                  if jb == 0 && kb == 0, continue; end  % Exclude the edge cell itself
                  mat.ijk_in_mixcell(mymx, :) = [i, nbdry - jb, nbdry - kb];
                  mymx = mymx + 1;
                  nghost_ylzl = nghost_ylzl + 1;
                end
              end

              % Handle yrzl edge
            elseif j == ncell_bdry(2) && k == nbdry
              for kb = 0:nbdry
                for jb = 0:nbdry
                  if jb == 0 && kb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [i, ncell_bdry(2) + jb, nbdry - kb];
                  mymx = mymx + 1;
                  nghost_yrzl = nghost_yrzl + 1;
                end
              end

              % Handle ylzr edge
            elseif j == nbdry && k == ncell_bdry(3)
              for kb = 0:nbdry
                for jb = 0:nbdry
                  if jb == 0 && kb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [i, nbdry - jb, ncell_bdry(3) + kb];
                  mymx = mymx + 1;
                  nghost_ylzr = nghost_ylzr + 1;
                end
              end

              % Handle xrzr edge
            elseif j == ncell_bdry(2) && k == ncell_bdry(3)
              for kb = 0:nbdry
                for jb = 0:nbdry
                  if jb == 0 && kb == 0, continue; end
                  mat.ijk_in_mixcell(mymx, :) = [i, ncell_bdry(2) + jb, ncell_bdry(3) + kb];
                  mymx = mymx + 1;
                  nghost_yrzr = nghost_yrzr + 1;
                end
              end
            end

            % Copy material properties for edge cells
            for idx = 1:nmix
              mymx = offset + idx;
              for m = 1:nm
                mat.matids_in_mixcell{mymx}(m) = mat.matids_in_mixcell{mx}(m);
                mat.vf_in_mixcell{mymx}(m) = mat.vf_in_mixcell{mx}(m);
                mat.rho_in_mixcell{mymx}(m) = mat.rho_in_mixcell{mx}(m);
                mat.pres_in_mixcell{mymx}(m) = mat.pres_in_mixcell{mx}(m);
                mat.ei_in_mixcell{mymx}(m) = mat.ei_in_mixcell{mx}(m);
              end
            end

            % Update the offset for the next edge
            offset = offset + nmix;

          elseif i == nbdry || i == ncell_bdry(1) || j == nbdry || j == ncell_bdry(2) || k == nbdry || k == ncell_bdry(3)  % face handling
            nmix = nbdry;
            nm = mat.nmat_in_mixcell(mx);
            mymx = offset;

            % Handle xl (left face)
            if i == nbdry
              for ib = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [nbdry - 1 - ib, j, k];
                mymx = mymx + 1;
              end

              % Handle xr (right face)
            elseif i == ncell_bdry(1)
              for ib = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [ncell_bdry(1) + 1 + ib, j, k];
                mymx = mymx + 1;
              end

              % Handle yl (bottom face)
            elseif j == nbdry
              for jb = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [i, nbdry - 1 - jb, k];
                mymx = mymx + 1;
                nghost_yl = nghost_yl + 1;
              end

              % Handle yr (top face)
            elseif j == ncell_bdry(2)
              for jb = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [i, ncell_bdry(2) + 1 + jb, k];
                mymx = mymx + 1;
                nghost_yr = nghost_yr + 1;
              end

              % Handle zl (front face)
            elseif k == nbdry
              for kb = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [i, j, nbdry - 1 - kb];
                mymx = mymx + 1;
                nghost_zl = nghost_zl + 1;
              end

              % Handle zr (back face)
            elseif k == ncell_bdry(3)
              for kb = 0:nbdry-1
                mat.ijk_in_mixcell(mymx, :) = [i, j, ncell_bdry(3) + 1 + kb];
                mymx = mymx + 1;
                nghost_zr = nghost_zr + 1;
              end
            end

            % Copy material properties for face cells
            for idx = 1:nmix
              mymx = offset + idx - 1;
              for m = 1:nm
                mat.matids_in_mixcell{mymx}(m) = mat.matids_in_mixcell{mx}(m);
                mat.vf_in_mixcell{mymx}(m) = mat.vf_in_mixcell{mx}(m);
                mat.rho_in_mixcell{mymx}(m) = mat.rho_in_mixcell{mx}(m);
                mat.pres_in_mixcell{mymx}(m) = mat.pres_in_mixcell{mx}(m);
                mat.ei_in_mixcell{mymx}(m) = mat.ei_in_mixcell{mx}(m);
              end
            end

            % Update the offset for the next face
            offset = offset + nmix;
          end

        end % for loop

      end % 2D or 3D

    end

    function bdry_cell_3d(~, mesh, btype_lower, btype_upper)
      % bdry_cell_3d Applies boundary conditions to 3D cell properties.
      %
      % This method handles transmitted boundary conditions for the cells on
      % the faces of a 3D grid. It updates material properties like number of
      % materials, material IDs, densities, internal energies, and pressures
      % for boundary cells.
      %
      % Inputs:
      %   mesh        - (c_mesh) The mesh object containing grid and material data.
      %   btype_lower - (array of int) Boundary type for the lower boundaries (xl, yl, zl).
      %   btype_upper - (array of int) Boundary type for the upper boundaries (xr, yr, zr).
      %
      % Uses:
      %   mesh.ncell_prob      - (array of int) Number of cells in each spatial dimension.
      %   mesh.nbdry_prob      - (int) Number of boundary cells.
      %   mesh.nmat_for_3dcell - (3D array of int) Number of materials in each 3D cell.
      %   mesh.matid_for_3dcell - (3D array of int) Material IDs for each 3D cell.
      %   mesh.rho_for_3dcell  - (3D array of double) Densities for each 3D cell.
      %   mesh.ei_for_3dcell   - (3D array of double) Internal energies for each 3D cell.
      %   mesh.pres_for_3dcell - (3D array of double) Pressures for each 3D cell.
      %
      % Modifies:
      %   mesh.nmat_for_3dcell - Updates the number of materials for the boundary cells.
      %   mesh.matid_for_3dcell - Updates the material IDs for the boundary cells.
      %   mesh.rho_for_3dcell  - Updates the densities for the boundary cells.
      %   mesh.ei_for_3dcell   - Updates the internal energies for the boundary cells.
      %   mesh.pres_for_3dcell - Updates the pressures for the boundary cells.
      %
      % Original C declaration:
      % void bdry_cell_3d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                   int ***nmat_for_3dcell, int ***matid_for_3dcell,
      %                   double ***rho_for_3dcell, double ***ei_for_3dcell, double ***pres_for_3dcell);

      ncell = mesh.ncell_prob;   % Number of cells in each direction
      nbdry = mesh.nbdry_prob;   % Number of boundary cells
      ncell_ext = ncell + 2 * nbdry;  % Extended cell size (including boundaries)
      nmat = mesh.nmat_mesh;       % Number of materials

      % xl edge without corners
      if btype_lower(1) == bdry_type.bdry_transmitted
        i0 = nbdry + 1;
        for k = nbdry + 1:ncell(3) + nbdry
          for j = nbdry + 1:ncell(2) + nbdry
            for i = 1:nbdry
              % Initialize cell properties
              mesh.rho_for_3dcell(k, j, i) = 0.0;
              mesh.ei_for_3dcell(k, j, i) = 0.0;
              mesh.pres_for_3dcell(k, j, i) = 0.0;

              % Loop over materials
              for m = 1:nmat
                % Update material properties
                mesh.vf_3dmat(k, j, i, m) = mesh.vf_3dmat(k, j, i0, m);
                mesh.rho_3dmat(k, j, i, m) = mesh.rho_3dmat(k, j, i0, m);
                mesh.ei_3dmat(k, j, i, m) = mesh.ei_3dmat(k, j, i0, m);
                mesh.pres_3dmat(k, j, i, m) = mesh.pres_3dmat(k, j, i0, m);

                % Aggregate properties for the cell
                mesh.rho_for_3dcell(k, j, i) = mesh.rho_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.rho_3dmat(k, j, i, m);
                mesh.ei_for_3dcell(k, j, i) = mesh.ei_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.ei_3dmat(k, j, i, m);
                mesh.pres_for_3dcell(k, j, i) = mesh.pres_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.pres_3dmat(k, j, i, m);
              end
            end
          end
        end
      end

      % xr face without edges
      if btype_lower(1) == bdry_type.bdry_transmitted
        i0 = ncell(1) + nbdry;
        for k = nbdry + 1:ncell(3) + nbdry
          for j = nbdry + 1:ncell(2) + nbdry
            for i = 1:nbdry
              i1 = i0 + i;
              % Initialize cell properties
              mesh.rho_for_3dcell(k, j, i1) = 0.0;
              mesh.ei_for_3dcell(k, j, i1) = 0.0;
              mesh.pres_for_3dcell(k, j, i1) = 0.0;

              % Loop over materials
              for m = 1:nmat
                % Update material properties
                mesh.vf_3dmat(k, j, i1, m) = mesh.vf_3dmat(k, j, i0, m);
                mesh.rho_3dmat(k, j, i1, m) = mesh.rho_3dmat(k, j, i0, m);
                mesh.ei_3dmat(k, j, i1, m) = mesh.ei_3dmat(k, j, i0, m);
                mesh.pres_3dmat(k, j, i1, m) = mesh.pres_3dmat(k, j, i0, m);

                % Aggregate properties for the cell
                mesh.rho_for_3dcell(k, j, i1) = mesh.rho_for_3dcell(k, j, i1) + ...
                  mesh.vf_3dmat(k, j, i1, m) * mesh.rho_3dmat(k, j, i1, m);
                mesh.ei_for_3dcell(k, j, i1) = mesh.ei_for_3dcell(k, j, i1) + ...
                  mesh.vf_3dmat(k, j, i1, m) * mesh.ei_3dmat(k, j, i1, m);
                mesh.pres_for_3dcell(k, j, i1) = mesh.pres_for_3dcell(k, j, i1) + ...
                  mesh.vf_3dmat(k, j, i1, m) * mesh.pres_3dmat(k, j, i1, m);
              end
            end
          end
        end
      end

      % yl face with edges
      if btype_lower(2) == bdry_type.bdry_transmitted
        j0 = nbdry + 1;
        for k = nbdry + 1:ncell(3) + nbdry
          for j = 1:nbdry
            for i = 1:ncell_ext(1)
              % Initialize cell properties
              mesh.rho_for_3dcell(k, j, i) = 0.0;
              mesh.ei_for_3dcell(k, j, i) = 0.0;
              mesh.pres_for_3dcell(k, j, i) = 0.0;

              % Loop over materials
              for m = 1:nmat
                % Update material properties
                mesh.vf_3dmat(k, j, i, m) = mesh.vf_3dmat(k, j0, i, m);
                mesh.rho_3dmat(k, j, i, m) = mesh.rho_3dmat(k, j0, i, m);
                mesh.ei_3dmat(k, j, i, m) = mesh.ei_3dmat(k, j0, i, m);
                mesh.pres_3dmat(k, j, i, m) = mesh.pres_3dmat(k, j0, i, m);

                % Aggregate properties for the cell
                mesh.rho_for_3dcell(k, j, i) = mesh.rho_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.rho_3dmat(k, j, i, m);
                mesh.ei_for_3dcell(k, j, i) = mesh.ei_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.ei_3dmat(k, j, i, m);
                mesh.pres_for_3dcell(k, j, i) = mesh.pres_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.pres_3dmat(k, j, i, m);
              end
            end
          end
        end
      end

      % yr face with edges
      if btype_lower(2) == bdry_type.bdry_transmitted
        j0 = ncell(2) + nbdry;
        for k = nbdry + 1:ncell(3) + nbdry
          for j = 1:nbdry
            j1 = j0 + j;
            for i = 1:ncell_ext(1)
              % Initialize cell properties
              mesh.rho_for_3dcell(k, j1, i) = 0.0;
              mesh.ei_for_3dcell(k, j1, i) = 0.0;
              mesh.pres_for_3dcell(k, j1, i) = 0.0;

              % Loop over materials
              for m = 1:nmat
                % Update material properties
                mesh.vf_3dmat(k, j1, i, m) = mesh.vf_3dmat(k, j0, i, m);
                mesh.rho_3dmat(k, j1, i, m) = mesh.rho_3dmat(k, j0, i, m);
                mesh.ei_3dmat(k, j1, i, m) = mesh.ei_3dmat(k, j0, i, m);
                mesh.pres_3dmat(k, j1, i, m) = mesh.pres_3dmat(k, j0, i, m);

                % Aggregate properties for the cell
                mesh.rho_for_3dcell(k, j1, i) = mesh.rho_for_3dcell(k, j1, i) + ...
                  mesh.vf_3dmat(k, j1, i, m) * mesh.rho_3dmat(k, j1, i, m);
                mesh.ei_for_3dcell(k, j1, i) = mesh.ei_for_3dcell(k, j1, i) + ...
                  mesh.vf_3dmat(k, j1, i, m) * mesh.ei_3dmat(k, j1, i, m);
                mesh.pres_for_3dcell(k, j1, i) = mesh.pres_for_3dcell(k, j1, i) + ...
                  mesh.vf_3dmat(k, j1, i, m) * mesh.pres_3dmat(k, j1, i, m);
              end
            end
          end
        end
      end

      % zl face including edges and corners
      if btype_lower(3) == bdry_type.bdry_transmitted
        k0 = nbdry + 1;
        for k = 1:nbdry
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              % Initialize cell properties
              mesh.rho_for_3dcell(k, j, i) = 0.0;
              mesh.ei_for_3dcell(k, j, i) = 0.0;
              mesh.pres_for_3dcell(k, j, i) = 0.0;

              % Loop over materials
              for m = 1:nmat
                % Update material properties
                mesh.vf_3dmat(k, j, i, m) = mesh.vf_3dmat(k0, j, i, m);
                mesh.rho_3dmat(k, j, i, m) = mesh.rho_3dmat(k0, j, i, m);
                mesh.ei_3dmat(k, j, i, m) = mesh.ei_3dmat(k0, j, i, m);
                mesh.pres_3dmat(k, j, i, m) = mesh.pres_3dmat(k0, j, i, m);

                % Aggregate properties for the cell
                mesh.rho_for_3dcell(k, j, i) = mesh.rho_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.rho_3dmat(k, j, i, m);
                mesh.ei_for_3dcell(k, j, i) = mesh.ei_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.ei_3dmat(k, j, i, m);
                mesh.pres_for_3dcell(k, j, i) = mesh.pres_for_3dcell(k, j, i) + ...
                  mesh.vf_3dmat(k, j, i, m) * mesh.pres_3dmat(k, j, i, m);
              end
            end
          end
        end
      end

      % zr face including edges and corners
      if btype_lower(3) == bdry_type.bdry_transmitted
        k0 = ncell(3) + nbdry;
        for k = 1:nbdry
          k1 = k0 + k;
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              % Initialize cell properties
              mesh.rho_for_3dcell(k1, j, i) = 0.0;
              mesh.ei_for_3dcell(k1, j, i) = 0.0;
              mesh.pres_for_3dcell(k1, j, i) = 0.0;

              % Loop over materials
              for m = 1:nmat
                % Update material properties
                mesh.vf_3dmat(k1, j, i, m) = mesh.vf_3dmat(k0, j, i, m);
                mesh.rho_3dmat(k1, j, i, m) = mesh.rho_3dmat(k0, j, i, m);
                mesh.ei_3dmat(k1, j, i, m) = mesh.ei_3dmat(k0, j, i, m);
                mesh.pres_3dmat(k1, j, i, m) = mesh.pres_3dmat(k0, j, i, m);

                % Aggregate properties for the cell
                mesh.rho_for_3dcell(k1, j, i) = mesh.rho_for_3dcell(k1, j, i) + ...
                  mesh.vf_3dmat(k1, j, i, m) * mesh.rho_3dmat(k1, j, i, m);
                mesh.ei_for_3dcell(k1, j, i) = mesh.ei_for_3dcell(k1, j, i) + ...
                  mesh.vf_3dmat(k1, j, i, m) * mesh.ei_3dmat(k1, j, i, m);
                mesh.pres_for_3dcell(k1, j, i) = mesh.pres_for_3dcell(k1, j, i) + ...
                  mesh.vf_3dmat(k1, j, i, m) * mesh.pres_3dmat(k1, j, i, m);
              end
            end
          end
        end
      end
    end

    function bdry_node_3d(~, mesh, btype_lower, btype_upper)
      % bdry_node_3d Applies boundary conditions to 3D node velocities.
      %
      % This method handles transmitted boundary conditions for node velocities
      % on the faces of a 3D grid, updating the velocity arrays for boundary
      % nodes, including the edges and corners.
      %
      % Inputs:
      %   mesh        - (c_mesh) The mesh object containing grid and velocity data.
      %   btype_lower - (array of int) Boundary type for the lower boundaries (xl, yl, zl).
      %   btype_upper - (array of int) Boundary type for the upper boundaries (xr, yr, zr).
      %
      % Uses:
      %   mesh.ncell_prob       - (array of int) Number of cells in each spatial dimension.
      %   mesh.nbdry_prob       - (int) Number of boundary cells.
      %   mesh.vel_for_3dnode   - (4D array) Velocities at each 3D node.
      %   mesh.vav_for_3dnode - (4D array) Averaged velocities at each 3D node.
      %
      % Modifies:
      %   mesh.vel_for_3dnode   - Updates the velocities at the boundary nodes.
      %   mesh.vav_for_3dnode - Updates the averaged velocities at the boundary nodes.
      %
      % Original C declaration:
      % void bdry_node_3d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                   double ****vel_for_3dnode, double ****vav_for_3dnode);

      dim = 3;
      ncell = mesh.ncell_prob;   % Number of cells in each direction
      nbdry = mesh.nbdry_prob;   % Number of boundary cells
      ncell_ext = ncell + 2 * nbdry;  % Extended cell size (including boundaries)
      nnode = ncell + 1;
      nnode_ext = ncell_ext + 1;
      
      % xl face without edges
      if btype_lower(1) == bdry_type.bdry_transmitted
        n0 = nbdry + 1;
        for k = nbdry + 1:nnode(3) + nbdry
          for j = nbdry + 1:nnode(2) + nbdry
            for i = 1:nbdry
              for dir = 1:dim
                mesh.vel_for_3dnode(k, j, i, dir) = mesh.vel_for_3dnode(k, j, n0, dir);
                mesh.vav_for_3dnode(k, j, i, dir) = mesh.vav_for_3dnode(k, j, n0, dir);
              end
            end
          end
        end
      end

      % xr face without edges
      if btype_upper(1) == bdry_type.bdry_transmitted
        n0 = nnode(1) + nbdry;
        for k = nbdry + 1:nnode(3) + nbdry
          for j = nbdry + 1:nnode(2) + nbdry
            for i = 1:nbdry
              n1 = n0 + i;
              for dir = 1:dim
                mesh.vel_for_3dnode(k, j, n1, dir) = mesh.vel_for_3dnode(k, j, n0, dir);
                mesh.vav_for_3dnode(k, j, n1, dir) = mesh.vav_for_3dnode(k, j, n0, dir);
              end
            end
          end
        end
      end

      % yl face with edges
      if btype_lower(2) == bdry_type.bdry_transmitted
        n0 = nbdry + 1;
        for k = nbdry + 1:nnode(3) + nbdry
          for j = 1:nbdry
            for i = 1:nnode_ext(1)
              for dir = 1:dim
                mesh.vel_for_3dnode(k, j, i, dir) = mesh.vel_for_3dnode(k, n0, i, dir);
                mesh.vav_for_3dnode(k, j, i, dir) = mesh.vav_for_3dnode(k, n0, i, dir);
              end
            end
          end
        end
      end

      % yr face with edges
      if btype_upper(2) == bdry_type.bdry_transmitted
        n0 = nnode(2) + nbdry;
        for k = nbdry + 1:nnode(3) + nbdry
          for j = 1:nbdry
            n1 = n0 + j;
            for i = 1:nnode_ext(1)
              for dir = 1:dim
                mesh.vel_for_3dnode(k, n1, i, dir) = mesh.vel_for_3dnode(k, n0, i, dir);
                mesh.vav_for_3dnode(k, n1, i, dir) = mesh.vav_for_3dnode(k, n0, i, dir);
              end
            end
          end
        end
      end

      % zl face including edges and corners
      if btype_lower(3) == bdry_type.bdry_transmitted
        n0 = nbdry + 1;
        for k = 1:nbdry
          for j = 1:nnode_ext(2)
            for i = 1:nnode_ext(1)
              for dir = 1:dim
                mesh.vel_for_3dnode(k, j, i, dir) = mesh.vel_for_3dnode(n0, j, i, dir);
                mesh.vav_for_3dnode(k, j, i, dir) = mesh.vav_for_3dnode(n0, j, i, dir);
              end
            end
          end
        end
      end

      % zr face including edges and corners
      if btype_upper(3) == bdry_type.bdry_transmitted
        n0 = nnode(3) + nbdry;
        for k = 1:nbdry
          n1 = n0 + k;
          for j = 1:nnode_ext(2)
            for i = 1:nnode_ext(1)
              for dir = 1:dim
                mesh.vel_for_3dnode(n1, j, i, dir) = mesh.vel_for_3dnode(n0, j, i, dir);
                mesh.vav_for_3dnode(n1, j, i, dir) = mesh.vav_for_3dnode(n0, j, i, dir);
              end
            end
          end
        end
      end
    end

    function bdry_cell_vel_2d(~, mesh, btype_lower, btype_upper)
      % bdry_cell_vel_2d Applies boundary conditions to 2D cell velocities.
      %
      % This function handles transmitted boundary conditions for the velocities
      % on the edges of a 2D grid. It modifies the velocity arrays for boundary
      % cells.
      %
      % Inputs:
      %   obj         - (bdry) The boundary object instance.
      %   mesh        - (c_mesh) The mesh object instance.
      %   btype_lower - (array) Boundary type for the lower boundary (xl, yl).
      %   btype_upper - (array) Boundary type for the upper boundary (xr, yr).

      dim = 2;
      ncell = mesh.ncell_prob;    % Number of cells in each spatial dimension
      nbdry = mesh.nbdry_prob;    % Number of boundary cells
      ncell_ext = ncell + 2 * nbdry;  % Number of cells including boundaries

      % Access mesh array for velocities
      vel_for_2dcell = mesh.vel_for_2dcell;

      % xl edge without corners
      if btype_lower(1) == bdry_type.bdry_transmitted
        i0 = nbdry + 1;
        for j = nbdry + 1:ncell(2) + nbdry
          for i = 1:nbdry
            for dir = 1:dim
              mesh.vel_for_2dcell(j, i, dir) = mesh.vel_for_2dcell(j, i0, dir);
            end
          end
        end
      end

      % xr edge without corners
      if btype_upper(1) == bdry_type.bdry_transmitted
        i0 = ncell(1) + nbdry;
        for j = nbdry + 1:ncell(2) + nbdry
          for i = 1:nbdry
            i1 = i0 + i;
            for dir = 1:dim
              mesh.vel_for_2dcell(j, i1, dir) = mesh.vel_for_2dcell(j, i0, dir);
            end
          end
        end
      end

      % yl edge including corners
      if btype_lower(2) == bdry_type.bdry_transmitted
        j0 = nbdry + 1;
        for j = 1:nbdry
          for i = 1:ncell_ext(1)
            for dir = 1:dim
              mesh.vel_for_2dcell(j, i, dir) = mesh.vel_for_2dcell(j0, i, dir);
            end
          end
        end
      end

      % yr edge including corners
      if btype_upper(2) == bdry_type.bdry_transmitted
        j0 = ncell(2) + nbdry;
        for j = 1:nbdry
          j1 = j0 + j;
          for i = 1:ncell_ext(1)
            for dir = 1:dim
              mesh.vel_for_2dcell(j1, i, dir) = mesh.vel_for_2dcell(j0, i, dir);
            end
          end
        end
      end
    end

    function bdry_cell_vel_3d(~, mesh, btype_lower, btype_upper)
      % bdry_cell_vel_3d Applies boundary conditions to 3D cell velocities.
      %
      % This method handles transmitted boundary conditions for cell velocities
      % on the faces of a 3D grid. It updates the velocity arrays for boundary
      % cells, including the edges and corners.
      %
      % Inputs:
      %   mesh        - (c_mesh) The mesh object containing grid and velocity data.
      %   btype_lower - (array of int) Boundary type for the lower boundaries (xl, yl, zl).
      %   btype_upper - (array of int) Boundary type for the upper boundaries (xr, yr, zr).
      %
      % Uses:
      %   mesh.ncell_prob      - (array of int) Number of cells in each spatial dimension.
      %   mesh.nbdry_prob      - (int) Number of boundary cells.
      %   mesh.vel_for_3dcell  - (4D array) Velocities at each 3D cell.
      %
      % Modifies:
      %   mesh.vel_for_3dcell  - Updates the velocities at the boundary cells.
      %
      % Original C declaration:
      % void bdry_cell_vel_3d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                       double ****vel_for_3dcell);

      dim = 3;
      ncell = mesh.ncell_prob;   % Number of cells in each direction
      nbdry = mesh.nbdry_prob;   % Number of boundary cells
      ncell_ext = ncell + 2 * nbdry;  % Extended cell size (including boundaries)

      % xl face without edges
      if btype_lower(1) == bdry_type.bdry_transmitted
        i0 = nbdry + 1;
        for k = nbdry+1:ncell(3)+nbdry
          for j = nbdry+1:ncell(2)+nbdry
            for i = 1:nbdry
              for dir = 1:dim
                mesh.vel_for_3dcell(k, j, i, dir) = mesh.vel_for_3dcell(k, j, i0, dir);
              end
            end
          end
        end
      end

      % xr face without edges
      if btype_upper(1) == bdry_type.bdry_transmitted
        i0 = ncell(1) + nbdry;
        for k = nbdry+1:ncell(3)+nbdry
          for j = nbdry+1:ncell(2)+nbdry
            for i = 1:nbdry
              i1 = i0 + i;
              for dir = 1:dim
                mesh.vel_for_3dcell(k, j, i1, dir) = mesh.vel_for_3dcell(k, j, i0, dir);
              end
            end
          end
        end
      end

      % yl face with edges
      if btype_lower(2) == bdry_type.bdry_transmitted
        j0 = nbdry + 1;
        for k = nbdry+1:ncell(3)+nbdry
          for j = 1:nbdry
            for i = 1:ncell_ext(1)
              for dir = 1:dim
                mesh.vel_for_3dcell(k, j, i, dir) = mesh.vel_for_3dcell(k, j0, i, dir);
              end
            end
          end
        end
      end

      % yr face with edges
      if btype_upper(2) == bdry_type.bdry_transmitted
        j0 = ncell(2) + nbdry;
        for k = nbdry+1:ncell(3)+nbdry
          for j = 1:nbdry
            j1 = j0 + j;
            for i = 1:ncell_ext(1)
              for dir = 1:dim
                mesh.vel_for_3dcell(k, j1, i, dir) = mesh.vel_for_3dcell(k, j0, i, dir);
              end
            end
          end
        end
      end

      % zl face including edges and corners
      if btype_lower(3) == bdry_type.bdry_transmitted
        k0 = nbdry + 1;
        for k = 1:nbdry
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              for dir = 1:dim
                mesh.vel_for_3dcell(k, j, i, dir) = mesh.vel_for_3dcell(k0, j, i, dir);
              end
            end
          end
        end
      end

      % zr face including edges and corners
      if btype_upper(3) == bdry_type.bdry_transmitted
        k0 = ncell(3) + nbdry;
        for k = 1:nbdry
          k1 = k0 + k;
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              for dir = 1:dim
                mesh.vel_for_3dcell(k1, j, i, dir) = mesh.vel_for_3dcell(k0, j, i, dir);
              end
            end
          end
        end
      end
    end

    function bdry_cell_1var_3d(obj, mesh, btype_lower, btype_upper)
      % bdry_cell_1var_3d Applies boundary conditions for a 3D cell array of materials.
      %
      % This function applies transmitted boundary conditions for the number of
      % materials in a 3D mesh, handling the lower and upper boundaries for
      % the x, y, and z directions, with specific rules for faces and edges.
      %
      % Inputs:
      %   mesh          - (c_mesh) The mesh object containing the number of cells and boundary data.
      %   btype_lower   - (array of int) Boundary type for the lower boundaries (xl, yl, zl).
      %   btype_upper   - (array of int) Boundary type for the upper boundaries (xr, yr, zr).
      %
      % Outputs:
      %   None. The function modifies the number of materials in each cell at the boundaries.
      %
      % Uses:
      %   mesh.ncell_prob      - (array of int) Number of cells in each direction [size: 1x3].
      %   mesh.nbdry_prob      - (int) Number of boundary layers.
      %
      % Modifies:
      %   obj.nmat_for_cell - (3D array of int) Updates the number of materials in cells
      %                          at the boundaries (xl, xr, yl, yr, zl, zr).
      %
      % Original C declaration:
      % void bdry_cell_1var_3d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                        int ***nmat_for_3dcell)

      dim = 3;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      ncell_ext = ncell + 2 * nbdry;  % Extended number of cells including boundary layers

      % xl face without edges
      if btype_lower(1) == bdry_type.bdry_transmitted
        i0 = nbdry + 1;
        for k = nbdry+1:ncell(3) + nbdry
          for j = nbdry+1:ncell(2) + nbdry
            for i = 1:nbdry
              obj.nmat_for_cell(k, j, i) = obj.nmat_for_cell(k, j, i0);
            end
          end
        end
      end

      % xr face without edges
      if btype_upper(1) == bdry_type.bdry_transmitted
        i0 = ncell(1) + nbdry;
        for k = nbdry+1:ncell(3) + nbdry
          for j = nbdry+1:ncell(2) + nbdry
            for i = 1:nbdry
              i1 = i0 + i;
              obj.nmat_for_cell(k, j, i1) = obj.nmat_for_cell(k, j, i0);
            end
          end
        end
      end

      % yl face with edges
      if btype_lower(2) == bdry_type.bdry_transmitted
        j0 = nbdry + 1;
        for k = nbdry+1:ncell(3) + nbdry
          for j = 1:nbdry
            for i = 1:ncell_ext(1)
              obj.nmat_for_cell(k, j, i) = obj.nmat_for_cell(k, j0, i);
            end
          end
        end
      end

      % yr face with edges
      if btype_upper(2) == bdry_type.bdry_transmitted
        j0 = ncell(2) + nbdry;
        for k = nbdry+1:ncell(3) + nbdry
          for j = 1:nbdry
            j1 = j0 + j;
            for i = 1:ncell_ext(1)
              obj.nmat_for_cell(k, j1, i) = obj.nmat_for_cell(k, j0, i);
            end
          end
        end
      end

      % zl face including edges and corners
      if btype_lower(3) == bdry_type.bdry_transmitted
        k0 = nbdry + 1;
        for k = 1:nbdry
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              obj.nmat_for_cell(k, j, i) = obj.nmat_for_cell(k0, j, i);
            end
          end
        end
      end

      % zr face including edges and corners
      if btype_upper(3) == bdry_type.bdry_transmitted
        k0 = ncell(3) + nbdry;
        for k = 1:nbdry
          k1 = k0 + k;
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              obj.nmat_for_cell(k1, j, i) = obj.nmat_for_cell(k0, j, i);
            end
          end
        end
      end
    end

    function bdry_cell_ragged_3d(obj, mesh, btype_lower, btype_upper)
      % bdry_cell_ragged_3d Applies boundary conditions for a 3D ragged cell array of materials.
      %
      % This method handles transmitted boundary conditions for the number of
      % materials and their properties (volume, mass, energy) in a 3D mesh,
      % updating boundary cells by copying data from inner cells.
      %
      % Inputs:
      %   obj           - (c_bdry) The current boundary object instance.
      %   mesh          - (c_mesh) The mesh object containing the cell data.
      %   btype_lower   - (array of int) Boundary type for the lower boundaries (xl, yl, zl).
      %   btype_upper   - (array of int) Boundary type for the upper boundaries (xr, yr, zr).
      %
      % Outputs:
      %   None. The function modifies the material properties at the boundaries.
      %
      % Uses:
      %   mesh.ncell_prob       - (array of int) Number of cells in each direction [size: 1x3].
      %   mesh.nbdry_prob       - (int) Number of boundary layers.
      %
      % Modifies:
      %   obj.nmat_for_cell     - (3D array of double) Number of materials for each cell.
      %   obj.matid_for_cell    - (cell array of int arrays) Material IDs for each cell.
      %   obj.vol_for_cell      - (cell array of double arrays) Volume for each material in each cell.
      %   obj.mass_for_cell     - (cell array of double arrays) Mass for each material in each cell.
      %   obj.ener_for_cell     - (cell array of double arrays) Energy for each material in each cell.
      %
      % Original C declaration:
      % void bdry_cell_ragged_3d(int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %                          int ***nmat_for_cell, int ****matid_for_cell, double ****vol_for_cell,
      %                          double ****mass_for_cell, double ****ener_for_cell)

      dim = 3;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      ncell_ext = ncell + 2 * nbdry;  % Extended number of cells including boundary layers

      % xl face without edges and corners
      if btype_lower(1) == bdry_type.bdry_transmitted
        i0 = nbdry + 1;
        for k = nbdry+1:ncell(3) + nbdry
          for j = nbdry+1:ncell(2) + nbdry
            for i = 1:nbdry
              obj.nmat_for_cell{k, j, i} = obj.nmat_for_cell{k, j, i0};
              nm = obj.nmat_for_cell{k, j, i0};
              for idx = 1:nm
                obj.matid_for_cell{k, j, i}{idx} = obj.matid_for_cell{k, j, i0}{idx};
                obj.vol_for_cell{k, j, i}{idx} = obj.vol_for_cell{k, j, i0}{idx};
                obj.mass_for_cell{k, j, i}{idx} = obj.mass_for_cell{k, j, i0}{idx};
                obj.ener_for_cell{k, j, i}{idx} = obj.ener_for_cell{k, j, i0}{idx};
              end
            end
          end
        end
      end

      % xr face without edges and corners
      if btype_upper(1) == bdry_type.bdry_transmitted
        i0 = ncell(1) + nbdry;
        for k = nbdry+1:ncell(3) + nbdry
          for j = nbdry+1:ncell(2) + nbdry
            for i = 1:nbdry
              i1 = i0 + i;
              obj.nmat_for_cell{k, j, i1} = obj.nmat_for_cell{k, j, i0};
              nm = obj.nmat_for_cell{k, j, i0};
              for idx = 1:nm
                obj.matid_for_cell{k, j, i1}{idx} = obj.matid_for_cell{k, j, i0}{idx};
                obj.vol_for_cell{k, j, i1}{idx} = obj.vol_for_cell{k, j, i0}{idx};
                obj.mass_for_cell{k, j, i1}{idx} = obj.mass_for_cell{k, j, i0}{idx};
                obj.ener_for_cell{k, j, i1}{idx} = obj.ener_for_cell{k, j, i0}{idx};
              end
            end
          end
        end
      end

      % yl face including edges
      if btype_lower(2) == bdry_type.bdry_transmitted
        j0 = nbdry + 1;
        for k = nbdry+1:ncell(3) + nbdry
          for j = 1:nbdry
            for i = 1:ncell_ext(1)
              obj.nmat_for_cell{k, j, i} = obj.nmat_for_cell{k, j0, i};
              nm = obj.nmat_for_cell{k, j0, i};
              for idx = 1:nm
                obj.matid_for_cell{k, j, i}{idx} = obj.matid_for_cell{k, j0, i}{idx};
                obj.vol_for_cell{k, j, i}{idx} = obj.vol_for_cell{k, j0, i}{idx};
                obj.mass_for_cell{k, j, i}{idx} = obj.mass_for_cell{k, j0, i}{idx};
                obj.ener_for_cell{k, j, i}{idx} = obj.ener_for_cell{k, j0, i}{idx};
              end
            end
          end
        end
      end

      % yr face including edges
      if btype_upper(2) == bdry_type.bdry_transmitted
        j0 = ncell(2) + nbdry;
        for k = nbdry+1:ncell(3) + nbdry
          for j = 1:nbdry
            j1 = j0 + j;
            for i = 1:ncell_ext(1)
              obj.nmat_for_cell{k, j1, i} = obj.nmat_for_cell{k, j0, i};
              nm = obj.nmat_for_cell{k, j0, i};
              for idx = 1:nm
                obj.matid_for_cell{k, j1, i}{idx} = obj.matid_for_cell{k, j0, i}{idx};
                obj.vol_for_cell{k, j1, i}{idx} = obj.vol_for_cell{k, j0, i}{idx};
                obj.mass_for_cell{k, j1, i}{idx} = obj.mass_for_cell{k, j0, i}{idx};
                obj.ener_for_cell{k, j1, i}{idx} = obj.ener_for_cell{k, j0, i}{idx};
              end
            end
          end
        end
      end

      % zl face including edges and corners
      if btype_lower(3) == bdry_type.bdry_transmitted
        k0 = nbdry + 1;
        for k = 1:nbdry
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              obj.nmat_for_cell{k, j, i} = obj.nmat_for_cell{k0, j, i};
              nm = obj.nmat_for_cell{k0, j, i};
              for idx = 1:nm
                obj.matid_for_cell{k, j, i}{idx} = obj.matid_for_cell{k0, j, i}{idx};
                obj.vol_for_cell{k, j, i}{idx} = obj.vol_for_cell{k0, j, i}{idx};
                obj.mass_for_cell{k, j, i}{idx} = obj.mass_for_cell{k0, j, i}{idx};
                obj.ener_for_cell{k, j, i}{idx} = obj.ener_for_cell{k0, j, i}{idx};
              end
            end
          end
        end
      end

      % zr face including edges and corners
      if btype_upper(3) == bdry_type.bdry_transmitted
        k0 = ncell(3) + nbdry;
        for k = 1:nbdry
          k1 = k0 + k;
          for j = 1:ncell_ext(2)
            for i = 1:ncell_ext(1)
              obj.nmat_for_cell{k1, j, i} = obj.nmat_for_cell{k0, j, i};
              nm = obj.nmat_for_cell{k0, j, i};
              for idx = 1:nm
                obj.matid_for_cell{k1, j, i}{idx} = obj.matid_for_cell{k0, j, i}{idx};
                obj.vol_for_cell{k1, j, i}{idx} = obj.vol_for_cell{k0, j, i}{idx};
                obj.mass_for_cell{k1, j, i}{idx} = obj.mass_for_cell{k0, j, i}{idx};
                obj.ener_for_cell{k1, j, i}{idx} = obj.ener_for_cell{k0, j, i}{idx};
              end
            end
          end
        end
      end
    end

  end
end
