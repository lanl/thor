classdef SMesh < handle

  properties
    % Problem definition variables
    dim_prob = 0;            % the dimensionality
    nbdry_prob = 0;          % the number of layers, typically 1
    ncell_prob = [0, 0, 0];  % the number of cells in three dimensions
    xl_prob = [0.0, 0.0, 0.0];  % the left ends of the simulation domain, excluding the ghost cells
    xr_prob = [0.0, 0.0, 0.0];  % the right ends of the simulation domain, excluding the ghost cells
    dx = [1.0, 1.0, 1.0];       % grid cell size

    nmat_mesh = 0;            % the number of materials in the mesh
    matids_mesh = [];         % material IDs in the mesh

    % Mesh variables (2D)
    rho_for_2dcell = [];   % density for each 2D cell
    ei_for_2dcell = [];    % internal energy density, per volume for each 2D cell
    pres_for_2dcell = [];  % pressure for each 2D cell
    vel_for_2dnode = [];   % velocity for each 2D node
    vav_for_2dnode = [];   % time-averaged velocity for advection for each 2D node

      % Mesh variables (3D)
      rho_for_3dcell = [];   % density for each 3D cell
      ei_for_3dcell = [];    % internal energy density, per volume for each 3D cell
      pres_for_3dcell = [];  % pressure for each 3D cell
      vel_for_3dnode = [];   % velocity for each 3D node
      vav_for_3dnode = [];   % time-averaged velocity for advection for each 3D node

    % Additional variables (2D) (Added by me - OK)
    vel_for_2dcell = [];   % velocity for each 2D *cell*
    cs_for_2dcell = [];    % sound speed for each 2D cell
    rho_for_2dcell_old = []; % density from the previous timestep
    es_for_2dcell_old = [];  % internal energy density from the previous timestep

      vel_for_3dcell = [];   % velocity for each 3D *cell*
      cs_for_3dcell = [];    % sound speed for each 3D cell
      rho_for_3dcell_old = []; % density from the previous timestep
      es_for_3dcell_old = [];  % internal energy density from the previous timestep

    % Material variables (2D)
    nmat_for_2dcell  = [];   % the number of materials in each 2D cell
    matid_for_2dcell = [];   % material ID if the cell is clean, otherwise -1

    % Material variables (3D)
    nmat_for_3dcell  = [];   % the number of materials in each 3D cell
    matid_for_3dcell = [];   % material ID if the cell is clean, otherwise -1

    % Mixed cell indices
    mixcell_for_2dcell = [];  % 2D cell to the index of mixed cells
    mixcell_for_3dcell = [];  % 2D cell to the index of mixed cells

    %extra
    vf_2dmat = [];
    rho_2dmat = [];
    ei_2dmat = [];
    pres_2dmat = [];

    vf_3dmat = [];
    rho_3dmat = [];
    ei_3dmat = [];
    pres_3dmat = [];

  end

  methods

    %%%%%%%%%%%%%%% implemented and tested
    function gradv = cal_cell_zgrad2d(~, dx, var)
      % cal_cell_zgrad2d - Calculate the gradient of a variable in a 2D cell.
      %
      % This function calculates the gradient of a variable across a 2D cell,
      % given the values of the variable at the cell's corners and the grid
      % spacing. The gradient is computed along the x and y directions.
      %
      % Inputs:
      %   dx    - (double array) The grid spacing in the x and y directions,
      %           size: [1, 2].
      %   var   - (double array) The values of the variable at the 3x3 grid points
      %           surrounding the cell, size: [3, 3].
      %
      % Outputs:
      %   gradv - (double array) The calculated gradient vector, with the x
      %           component in gradv(1) and the y component in gradv(2),
      %           size: [1, 2].
      %
      % Original C declaration:
      % void cal_cell_zgrad2d(double *dx, double var[3][3], double *gradv)

      % Initialize variables
      gradv = zeros(1, 2);
      v0 = var(2, 2);  % MATLAB uses 1-based indexing, so var[1][1] becomes var(2,2)

      % Calculate gradient in the x-direction
      dfsum = 0.0;
      for j = 1:3
        dfsum = dfsum + (var(j, 3) - v0) - (var(j, 1) - v0);
      end
      gradv(1) = dfsum / (6.0 * dx(1));

      % Calculate gradient in the y-direction
      dfsum = 0.0;
      for i = 1:3
        dfsum = dfsum + (var(3, i) - v0) - (var(1, i) - v0);
      end
      gradv(2) = dfsum / (6.0 * dx(2));
    end

    %%%%%%%%%%%%%%% implemented but not tested
    function [xl, xr, ncell, nbdry, nmat_ea_2dcell, matid_ea_2dcell,...
        mixcell_ea_2dcell, nmat_ea_3dcell, matid_ea_3dcell,...
        mixcell_ea_3dcell] = mesh_pass_mat(obj)

      xl = obj.xl_prob;
      xr = obj.xr_prob;
      ncell = obj.ncell_prob;
      nbdry = obj.nbdry_prob;
      nmat_ea_2dcell = obj.nmat_for_2dcell;
      matid_ea_2dcell = obj.matid_for_2dcell;
      mixcell_ea_2dcell = obj.mixcell_for_2dcell;

      nmat_ea_3dcell = obj.nmat_for_3dcell;
      matid_ea_3dcell = obj.matid_for_3dcell;
      mixcell_ea_3dcell = obj.mixcell_for_3dcell;

    end

    function courant_cs = courant_from_cs(obj, dt)
      % courant_from_cs Calculate the maximum Courant number based on sound speed.
      %
      % This method computes the maximum Courant number across all cells in a 2D
      % or 3D computational mesh. The Courant number is crucial for determining
      % the stability of numerical simulations involving time-dependent PDEs.
      %
      % Inputs:
      %   obj        - (c_mesh) Instance of the c_mesh class, containing
      %                information about the computational mesh.
      %   dt         - (double) The time step size used in the simulation.
      %
      % Outputs:
      %   courant_cs - (double) The maximum Courant number computed over all
      %                cells in the mesh.
      %
      % Uses:
      %   obj.dim_prob       - (int) Dimensionality of the problem (2 or 3).
      %   obj.dx             - (array of double) Grid spacing in each spatial
      %                        dimension.
      %   obj.cs_for_2dcell  - (2D array of double) Sound speed array for 2D
      %                        cells.
      %   obj.cs_for_3dcell  - (3D array of double) Sound speed array for 3D
      %                        cells.
      %
      % Original C declaration:
      %   void courant_from_cs(int dim, int *ncell, int nbdry, double *dx,
      %                        double **cs_for_2dcell, double ***cs_for_3dcell,
      %                        double dt, double *courant_cs);
      dim = obj.dim_prob;
      dxmin = min(obj.dx(1:dim));

      if dim == 2
        courant_cs = max(obj.cs_for_2dcell, [], "all");
      elseif dim == 3
        courant_cs = max(obj.cs_for_3dcell, [], "all");
      end
      courant_cs = courant_cs * dt / dxmin;
    end

    function courant_vel = courant_from_vel(obj, dt)
      % courant_from_vel Calculate the maximum Courant number based on sound speed.
      %
      % This method computes the maximum Courant number across all cells in a 2D
      % or 3D computational mesh. The Courant number is crucial for determining
      % the stability of numerical simulations involving time-dependent PDEs.
      %
      % Inputs:
      %   obj        - (c_mesh) Instance of the c_mesh class, containing
      %                information about the computational mesh.
      %   dt         - (double) The time step size used in the simulation.
      %
      % Outputs:
      %   courant_vel - (double) The maximum Courant number computed over all
      %                cells in the mesh.
      %
      % Uses:
      %   obj.dim_prob       - (int) Dimensionality of the problem (2 or 3).
      %   obj.dx             - (array of double) Grid spacing in each spatial
      %                        dimension.
      % Original C declaration:
      %   void courant_from_vel(int dim, int *ncell, int nbdry, double *dx,
      %                        double ***vel_2dnode, double ****vel_3dnode,
      %                        double dt, double *courant_cs);
      dim = obj.dim_prob;
      dx = obj.dx(1:dim);
      ncell_ext = obj.ncell_prob + 2*obj.nbdry_prob;

      courant_vel = 0.0;
      if dim == 2
          for j = 1:ncell_ext(2)
              for i = 1:ncell_ext(1)
                  for dir = 1:dim
                     c = abs(obj.vel_for_2dnode(j, i, dir))*dt/dx(dir);
                     courant_vel = max(courant_vel, c);
                  end
              end
          end
      elseif dim == 3
          for k = 1:ncell_ext(3)
            for j = 1:ncell_ext(2)
              for i = 1:ncell_ext(1)
                  for dir = 1:dim
                     c = abs(obj.vel_for_3dnode(k, j, i, dir))*dt/dx(dir);
                     courant_vel = max(courant_vel, c);
                  end
              end
            end
          end
      end
    end

    function set_mesh(obj, dim, xl, xr, ncell, nbdry, nmat, matids)
      % Function: set_mesh
      %
      % Description:
      %   Initializes and sets up the mesh for the simulation, including the
      %   domain boundaries, cell information, and mesh variables for both
      %   2D and 3D grids. Depending on the dimensionality (2D or 3D),
      %   it allocates memory for various mesh-related variables, such as
      %   the number of materials per cell, material IDs, mixed cell indices,
      %   and physical quantities like density, internal energy, pressure,
      %   and velocity.
      %
      % Inputs:
      %   obj   - (c_mesh) The instance of the c_mesh class.
      %   dim   - (int) Dimensionality of the problem (2 or 3).
      %   xl    - (double array) Lower bounds of the simulation domain in each dimension.
      %   xr    - (double array) Upper bounds of the simulation domain in each dimension.
      %   ncell - (int array) Number of cells in each dimension.
      %   nbdry - (int) Number of boundary layers (typically 1).
      %   nmat  - (int) Number of materials in the mesh.
      %   matids- (int array) Material IDs in the mesh.
      %
      % Outputs:
      %   None. (The function modifies the properties of the obj instance directly.)

      assert((dim > 1) && (dim <= 3));

      % Initialize variables
      obj.xl_prob = zeros(1, dim);
      obj.xr_prob = zeros(1, dim);
      obj.ncell_prob = zeros(1, dim);
      ncell_ext = zeros(1, dim);
      nnode_ext = zeros(1, dim);

      obj.dim_prob = dim;
      lsize = 1;
      lsize_node = 1;

      % Set up domain and cell information
      for i = 1:dim
        obj.xl_prob(i) = xl(i);
        obj.xr_prob(i) = xr(i);
        obj.ncell_prob(i) = ncell(i);
        obj.dx(i) = (xr(i) - xl(i))/ncell(i);

        ncell_ext(i) = ncell(i) + 2 * nbdry;
        nnode_ext(i) = ncell_ext(i) + 1;
        lsize = lsize * ncell_ext(i);
        lsize_node = lsize_node * nnode_ext(i);
      end

      obj.nbdry_prob = nbdry;
      obj.nmat_mesh = nmat;
      obj.matids_mesh = matids;

      if dim == 2
        % 2D Arrays Initialization
        obj.nmat_for_2dcell = zeros(ncell_ext(2), ncell_ext(1));
        obj.matid_for_2dcell = zeros(ncell_ext(2), ncell_ext(1));
        obj.mixcell_for_2dcell = -ones(ncell_ext(2), ncell_ext(1));

        obj.rho_for_2dcell = zeros(ncell_ext(2), ncell_ext(1));
        obj.ei_for_2dcell = zeros(ncell_ext(2), ncell_ext(1));
        obj.pres_for_2dcell = zeros(ncell_ext(2), ncell_ext(1));

        % Node velocity and time-averaged velocity arrays
        obj.vel_for_2dnode = zeros(nnode_ext(2), nnode_ext(1), dim);
        obj.vav_for_2dnode = zeros(nnode_ext(2), nnode_ext(1), dim);

      elseif dim == 3
        % 3D Arrays Initialization
        obj.nmat_for_3dcell = ones(ncell_ext(3), ncell_ext(2), ncell_ext(1));
        obj.matid_for_3dcell = ones(ncell_ext(3), ncell_ext(2), ncell_ext(1));
        obj.mixcell_for_3dcell = -ones(ncell_ext(3), ncell_ext(2), ncell_ext(1));

        obj.rho_for_3dcell = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));
        obj.ei_for_3dcell = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));
        obj.pres_for_3dcell = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));

        % Node velocity and time-averaged velocity arrays
        obj.vel_for_3dnode = zeros(nnode_ext(3), nnode_ext(2), nnode_ext(1), dim);
        obj.vav_for_3dnode = zeros(nnode_ext(3), nnode_ext(2), nnode_ext(1), dim);
      end

    end

    function nmix_cell = set_mesh_mat(obj, dim, btype_lower, btype_upper, ...
        nmat, matids_ea_mat, is_solid, gamma_ea_mat, ...
        nreg, reg2matids, reg_shape, rho_ea_reg, ...
        pres_ea_reg, ei_ea_reg, v_ea_reg)
      % Function: set_mesh_mat
      %
      % Description:
      %   This method initializes the material properties for a computational
      %   mesh and sets up the material distribution across the mesh cells.
      %   The method handles both 2D and 3D meshes, determining the number of
      %   materials in each cell, assigning material IDs, and setting properties
      %   such as density, internal energy, and pressure. The method also handles
      %   the boundary conditions for the material properties.
      %
      % Inputs:
      %   obj - Instance of the c_mesh class
      %   dim - Dimensionality of the problem (2 or 3) [double]
      %   btype_lower - Boundary type at the lower end of each dimension [int]
      %   btype_upper - Boundary type at the upper end of each dimension [int]
      %   nmat - Number of materials in the problem [uint8]
      %   matids_ea_mat - Array of material IDs for each material [uint8]
      %   is_solid - Array indicating if each material is solid [uint8]
      %   gamma_ea_mat - Array of adiabatic index values for each material [double]
      %   nreg - Number of regions in the mesh [int]
      %   reg2matids - Array mapping regions to material IDs [int]
      %   reg_shape - Array of region_shape objects defining the shape of each
      %               region [region_shape]
      %   rho_ea_reg - Array of density values for each region [double]
      %   pres_ea_reg - Array of pressure values for each region [double]
      %   ei_ea_reg - Array of internal energy density values for each region [double]
      %   v_ea_reg - Cell array of velocity arrays for each region [cell]
      %
      % Outputs:
      %   nmix_cell - Number of mixed cells in the mesh [int]
      %
      % Notes:
      %   - This method updates the properties of the mesh object in place,
      %     such as material distributions and cell properties, based on the
      %     input material and region properties.
      %   - The method also calculates the number of mixed cells and applies
      %     boundary conditions to the cell variables.
      %   - The function handles the different cases for 2D and 3D meshes separately.

      global mat io xio;
      lsize = 1;
      dx = obj.dx(1:dim);
      ncell_ext = zeros(1, dim);

      % Calculate cell sizes and initialize arrays based on dimensionality
      for i = 1:dim
        ncell_ext(i) = obj.ncell_prob(i) + 2*obj.nbdry_prob;
        lsize = lsize*ncell_ext(i);
      end

      % Allocate memory for velocity arrays
      obj.vel_for_2dcell = zeros(1,1,1);
      obj.vel_for_3dcell = zeros(1,1,1,1);
      if dim == 2
        sizes = [ncell_ext(2), ncell_ext(1), dim];
        obj.vel_for_2dcell = zeros(sizes);
      elseif dim == 3
        sizes = [ncell_ext(3), ncell_ext(2), ncell_ext(1), dim];
        obj.vel_for_3dcell = zeros(sizes);
      end

      % Instantiate c_mat object and call set_mat method
      nmix_cell = mat.set_mat(obj, btype_lower, btype_upper, ...
        is_solid, gamma_ea_mat, nreg, reg2matids, reg_shape, rho_ea_reg, ...
        pres_ea_reg, ei_ea_reg, v_ea_reg);

      % Get node velocities and assign them to the appropriate properties
      if dim == 2
        obj.get_node_vel_2d();
      elseif dim == 3
        obj.get_node_vel_3d();
      end

      % Generate material polyhedrons
      mat.get_mpoly(obj);

      % Write mesh and material data to file
      t = 0.0;
      ncycle = 0;
      fileid = io.create_file("file_initial", t, ncycle);
      meshid = io.write_mesh_mat(fileid, "mesh", obj);
      io.write_mesh_vars(fileid, meshid, obj);

      xio.create_new_output_files();
      xio.write_dump(obj, mat, ncycle, t);

      mat.mat_write_mat(fileid, dim, -1);
      io.close_file(fileid);

      % Free allocated memory
      if dim == 2
        obj.vel_for_2dcell = [];
      elseif dim == 3
        obj.vel_for_3dcell = [];
      end

    end

    function get_node_vel_2d(obj)
      % get_node_vel_2d Computes the velocity at each node in a 2D mesh.
      %
      % This method calculates the velocity at each node of a 2D computational
      % mesh by averaging the velocities of the surrounding cells, weighted by
      % the mass (density) of each neighboring cell, over the total mass around
      % each node.
      %
      % Modifies:
      %   obj.vel_for_2dnode     - (array) Storing the computed node velocities.
      %
      % Uses:
      %   obj.ncell_prob         - (array) Number of cells in each spatial dimension.
      %   obj.nbdry_prob         - (int) Number of boundary cells.
      %   obj.rho_for_2dcell     - (2D array) Density for each cell.
      %   obj.vel_for_2dcell     - (3D array) Velocity for each cell in both
      %                            directions (x, y).
      %
      % See also: get_node_vel_3d

      global tiny;
      dim = 2;
      nnode_last = obj.ncell_prob + obj.nbdry_prob;

      for j = obj.nbdry_prob:nnode_last(2)
        for i = obj.nbdry_prob:nnode_last(1)
          mass = 0.0;
          % Initialize node velocities to zero
          for dir = 1:dim
            obj.vel_for_2dnode(j + 1, i + 1, dir) = 0.0;
          end
          % Accumulate contributions from neighboring cells
          for jc = j-1:j
            for ic = i-1:i
              for dir = 1:dim
                obj.vel_for_2dnode(j + 1, i + 1, dir) = ...
                  obj.vel_for_2dnode(j + 1, i + 1, dir) + ...
                  (obj.rho_for_2dcell(jc + 1, ic + 1) * obj.vel_for_2dcell(jc + 1, ic + 1, dir));
              end
              mass = mass + obj.rho_for_2dcell(jc + 1, ic + 1);
            end
          end
          % Normalize the velocities by the total mass
          massinv = 1.0 / (mass + tiny);
          for dir = 1:dim
            obj.vel_for_2dnode(j + 1, i + 1, dir) = obj.vel_for_2dnode(j + 1, i + 1, dir) * massinv;
          end
        end
      end
    end

    function get_node_vel_3d(obj)
      % get_node_vel_3d Computes the velocity at each node in a 3D mesh.
      %
      % This method calculates the velocity at each node of a 3D computational
      % mesh by averaging the velocities of the surrounding cells, weighted by
      % the mass (density) of each neighboring cell, over the total mass around
      % each node.
      %
      % Modifies:
      %   obj.vel_for_3dnode      - (array) storing the computed node velocities.
      %
      % Uses:
      %   obj.ncell_prob         - (array) Number of cells in each spatial dimension.
      %   obj.nbdry_prob         - (int) Number of boundary cells.
      %   obj.rho_for_3dcell     - (3D array) Density for each cell.
      %   obj.vel_for_3dcell     - (4D array) Velocity for each cell in all three
      %                            directions (x, y, z).
      %   obj.vel_for_3dnode     - (4D array) Velocity at each node, calculated and
      %                            stored by this function.
      %
      % See also: get_node_vel_2d

      global tiny;
      dim = 3;
      nnode_last = obj.ncell_prob + obj.nbdry_prob;

      for k = obj.nbdry_prob:nnode_last(3)
        for j = obj.nbdry_prob:nnode_last(2)
          for i = obj.nbdry_prob:nnode_last(1)
            mass = 0.0;
            % Initialize node velocities to zero
            for dir = 1:dim
              obj.vel_for_3dnode(k + 1, j + 1, i + 1, dir) = 0.0;
            end
            % Accumulate contributions from neighboring cells
            for kc = k-1:k
              for jc = j-1:j
                for ic = i-1:i
                  for dir = 1:dim
                    obj.vel_for_3dnode(k + 1, j + 1, i + 1, dir) = ...
                      obj.vel_for_3dnode(k + 1, j + 1, i + 1, dir) + ...
                      (obj.rho_for_3dcell(kc + 1, jc + 1, ic + 1) * ...
                      obj.vel_for_3dcell(kc + 1, jc + 1, ic + 1, dir));
                  end
                  mass = mass + obj.rho_for_3dcell(kc + 1, jc + 1, ic + 1);
                end
              end
            end
            % Normalize the velocities by the total mass
            massinv = 1.0 / (mass + tiny);
            for dir = 1:dim
              obj.vel_for_3dnode(k + 1, j + 1, i + 1, dir) = ...
                obj.vel_for_3dnode(k + 1, j + 1, i + 1, dir) * massinv;
            end
          end
        end
      end
    end

    function gradv = cal_cell_zgrad3d(~, dx, var)
      % cal_cell_zgrad3d - Calculate the gradient of a variable in a 3D cell.
      %
      % This function calculates the gradient of a variable across a 3D cell,
      % given the values of the variable at the cell's neighboring points and
      % the grid spacing. The gradient is computed along the x, y, and z directions.
      %
      % Inputs:
      %   ~      - Unused object reference, indicating that this function does not
      %            use any properties of the `c_mesh` class instance.
      %   dx     - (double array) The grid spacing in the x, y, and z directions,
      %            size: [1, 3].
      %   var    - (double array) The values of the variable at the 3x3x3 grid points
      %            surrounding the cell, size: [3, 3, 3].
      %
      % Outputs:
      %   gradv  - (double array) The calculated gradient vector, with the x, y,
      %            and z components in gradv(1), gradv(2), and gradv(3),
      %            respectively, size: [1, 3].
      %
      % Uses:
      %   This function does not use any properties from the `c_mesh` class instance
      %   or other objects.
      %
      % Modifies:
      %   This function does not modify any properties or external variables; it only
      %   computes and returns the gradient vector `gradv`.
      %
      % Original C declaration:
      % void cal_cell_zgrad3d(double *dx, double var[3][3][3], double *gradv)

      % Initialize the output gradient vector
      gradv = zeros(1, 3);

      % Extract the central value
      v0 = var(2, 2, 2);

      % Calculate the gradient in the x-direction
      dfsum = 0.0;
      for k = 1:3
        for j = 1:3
          dfsum = dfsum + (var(k, j, 3) - v0) - (var(k, j, 1) - v0);
        end
      end
      gradv(1) = dfsum / (18.0 * dx(1));

      % Calculate the gradient in the y-direction
      dfsum = 0.0;
      for k = 1:3
        for i = 1:3
          dfsum = dfsum + (var(k, 3, i) - v0) - (var(k, 1, i) - v0);
        end
      end
      gradv(2) = dfsum / (18.0 * dx(2));

      % Calculate the gradient in the z-direction
      dfsum = 0.0;
      for j = 1:3
        for i = 1:3
          dfsum = dfsum + (var(3, j, i) - v0) - (var(1, j, i) - v0);
        end
      end
      gradv(3) = dfsum / (18.0 * dx(3));
    end

    function mesh_sound_speed(obj, nmat, matids, is_solid, gamma_ea_mat)
      % mesh_sound_speed Calculates the sound speed for each cell in the mesh.
      %
      % This method computes the sound speed for each 2D or 3D cell in the mesh
      % based on the materials present in the cells. It distinguishes between
      % solid and non-solid materials, applying different methods to compute
      % the sound speed.
      %
      % Inputs:
      %   obj          - (c_mesh) The current mesh object instance.
      %   nmat         - (int) Number of materials in the mesh.
      %   matids       - (array of int) Material IDs for each material.
      %   is_solid     - (array of int) Indicator array for solid materials.
      %   gamma_ea_mat - (array of double) Gamma values for each material.
      %
      % Outputs:
      %   None. The function updates the sound speed arrays for 2D and 3D cells.
      %
      % Uses:
      %   obj.dim_prob          - (int) Dimensionality of the mesh (2D or 3D).
      %   obj.ncell_prob        - (array of int) Number of cells in each spatial dimension.
      %   obj.nbdry_prob        - (int) Number of boundary cells.
      %   obj.nmat_for_2dcell   - (2D array of int) Number of materials in each 2D cell.
      %   obj.matid_for_2dcell  - (2D array of int) Material IDs for each 2D cell.
      %   obj.rho_for_2dcell    - (2D array of double) Density for each 2D cell.
      %   obj.pres_for_2dcell   - (2D array of double) Pressure for each 2D cell.
      %   obj.nmat_for_3dcell   - (3D array of int) Number of materials in each 3D cell.
      %   obj.matid_for_3dcell  - (3D array of int) Material IDs for each 3D cell.
      %   obj.rho_for_3dcell    - (3D array of double) Density for each 3D cell.
      %   obj.pres_for_3dcell   - (3D array of double) Pressure for each 3D cell.
      %
      % Modifies:
      %   obj.cs_for_2dcell     - (2D array of double) Updates sound speed in each 2D cell.
      %   obj.cs_for_3dcell     - (3D array of double) Updates sound speed in each 3D cell.
      %
      % Original C declaration:
      % void mesh_sound_speed(int nmat, int *matids, int *is_solid, double *gamma_ea_mat,
      %              double **cs_for_2dcell,  double ***cs_for_3dcell);

      global eos tiny;  % Global object for EOS (Equation of State)

      dim = obj.dim_prob;       % Dimensionality of the mesh
      ncell = obj.ncell_prob;   % Number of cells in each direction
      nbdry = obj.nbdry_prob;   % Number of boundary cells

      % Initialize extended cell size
      ncell_ext = ncell + 2 * nbdry;

      % 2D case
      if dim == 2
          for j = 1:ncell_ext(2)
              for i = 1:ncell_ext(1)
                  obj.cs_for_2dcell(j,i) = 0;
                  for m = 1: nmat
                      if obj.vf_2dmat(j,i,m)>0 
                          if ~is_solid(m)
                              cs = sqrt(gamma_ea_mat(m) * obj.pres_2dmat(j,i,m) / ...
                                  (obj.rho_2dmat(j,i,m) + tiny));
                          else
                              cs = eos.sound_speed_solid(obj.rho_2dmat(j,i,m));
                          end
                          obj.cs_for_2dcell(j,i) = obj.cs_for_2dcell(j,i) + (obj.vf_2dmat(j,i,m)*cs);
                      end
                  end
              end
          end

      % 3D case
      elseif dim == 3
          for k = 1:ncell_ext(3)
              for j = 1:ncell_ext(2)
                  for i = 1:ncell_ext(1)
                      obj.cs_for_3dcell(k,j,i) = 0;
                      for m = 1: nmat
                          if obj.vf_3dmat(k,j,i,m)>0 
                              if ~is_solid(m)
                                  cs = sqrt(gamma_ea_mat(m) * obj.pres_3dmat(k,j,i,m) / ...
                                      (obj.rho_3dmat(k,j,i,m) + tiny));
                              else
                                  cs = eos.sound_speed_solid(obj.rho_3dmat(k,j,i,m));
                              end
                              obj.cs_for_3dcell(k,j,i) = obj.cs_for_3dcell(k,j,i) + (obj.vf_3dmat(k,j,i,m)*cs);
                          end
                      end
                  end
              end
          end
      end
    end

%%%%%%%%%%%%%%% empty stubs
    function cal_zgrad(obj, dim, var, zgrad)
      % TODO: Implement cal_zgrad
    end

    function cal_zgrad2d(obj, var, zgrad)
      % TODO: Implement cal_zgrad2d
    end

    function cal_zgrad3d(obj, var, zgrad)
      % TODO: Implement cal_zgrad3d
    end

    function mesh_pass_mesh_data(obj, nmat_ea_2dcell, ...
        matid_ea_2dcell, rho_ea_2dcell, pres_ea_2dcell, ...
        ei_ea_2dcell, vel_ea_2dnode, velav_ea_2dnode, ...
        nmat_ea_3dcell, matid_ea_3dcell, rho_ea_3dcell, ...
        pres_ea_3dcell, ei_ea_3dcell, vel_ea_3dnode, ...
        velav_ea_3dnode)
      % TODO: Implement mesh_pass_mesh_data
    end

    function mesh_pass_viz_file_ids(obj, fid, mid)
      % TODO: Implement mesh_pass_viz_file_ids
    end

  end
end

