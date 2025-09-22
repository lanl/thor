classdef c_mat < handle
  % c_mat Class for managing material properties in a computational mesh.
  %
  % This class is responsible for handling material properties in a computational
  % mesh. It stores information about materials, such as their densities,
  % pressures, internal energies, volume fractions, and velocities. The class
  % also includes methods for computing material properties for a specific cell
  % in the mesh, based on input regions and materials.
  %
  properties
    vfmin = 1.0e-06;        % smallest volume fraction, otherwise, considered zero.

    nmat_prob = 0;          % the number of materials in the problem.
    matids_prob = [];       % matids of each material

    nmixcell = 0;           % total number of mixed cells in the mesh
    nmixcell_int = 0;       % number of mixed cells in interior cells of the mesh
    nmat_in_mixcell = [];   % [1:nmixcell] number of materials in each mixed cell
    % see also: mesh.nmat_for_2dcell / _3dcell
    matids_in_mixcell = {}; % {1:nmixcell}[1:nm] material IDs of each material in each mixed cell
    % see also: mesh.matid_for_2dcell / _3dcell
    ijk_in_mixcell = [];    % [1:nmixcell,1:dim] mapping from the list of mixed cells to the mesh
    % the inverse mapping is stored in mesh.mixcell_for_2dcell / _3dcell

    vf_in_mixcell = {};     % {1:nmixcell}[1:nm] volume fractions of each material in each mixed cell
    rho_in_mixcell = {};    % {1:nmixcell}[1:nm] densities of each material in each mixed cell
    pres_in_mixcell = {};   % {1:nmixcell}[1:nm] pressures of each material in each mixed cell
    ei_in_mixcell = {};     % {1:nmixcell}[1:nm] internal energy density in each mixed cell
    vfgrad_mixcell = {};    % {1:nmixcell}[1:nm,1:dim] volume fraction gradient for e.m. in e.m.c.

    nmixcell_mpoly = 0;     % number of mixed cells in which there are material polygons
    mix_mpoly_to_mix = [];  % [1:nmixcell_mpoly] mapping from shortened index of mixed cells
    % with material polygons into the full index:
    % [1:nmixcell_mpoly] -> [1,nmixcell]
    mix_to_mpoly_mix = [];  % inverse mapping: [1,nmixcell] -> [1,nmixcell_mpoly]

    nnode_in_mixcell = [];  % [1:nmixcell] total number of material polygon nodes in each mixcell
    coords_in_mixcell = {}; % {1:nmixcell}[1:nnodes,1:dim] coordinates of nodes in each mixcell

    nnode_for_minterface_in_mixcell = {}; % {1:nmixcell}[1:nm-1] - how many nodes are there on each
    % of the material interfaces (in 2D, always 2)
    % Note: there are exactly nm - 1 interfaces for nm materials
    nodes_for_minterface_in_mixcell = {}; % {1:nmixcell}{1:nm-1}[1:nn] - list of nodes for every
    % interface in the mixcell. Coordinates of the nodes are
    % stored in `coords_in_mixcell`
    % Note: nested cell array!

    % For 2D
    nnode_for_mpoly_in_mixcell = {}; % {1:nmixcell}[whichpolygon:1..nm] - number of nodes in each material polygon
    % number of material polygons is the same as the number of materials, i.e.
    % nmat_in_mixcell[mx]
    nodes_for_mpoly_in_mixcell = {}; % {1:nmixcell}{whichpolygon}[nodelist] - list of nodes for each material
    % polygon. Nodes are specified in nnode_in_mixcell and coords_in_mixcell

    % For 3D
    nface_for_mpoly_in_mixcell = {};
    nnode_for_face_ea_mpoly_in_mixcell = {};
    nodelist_for_face_ea_mpoly_in_mixcell = {};

  end

  methods

    %%%%%%%%%%%%%%% implemented but not tested
    function set_mpoly_storage(obj, dim)
      % set_mpoly_storage - Allocate or clear storage for material polygons and
      %                     interfaces based on the dimensionality and number of
      %                     mixed cells.
      %
      % This function initializes or clears the storage for material polygons
      % (2D/3D) and interfaces in mixed cells. Depending on the dimensionality,
      % it allocates storage for 2D polygons or 3D polyhedrons and interfaces.
      %
      % Inputs:
      %   obj     - (c_mat object) The instance of the c_mat class containing
      %             material data and mixed cell information.
      %   dim     - (int) The dimensionality of the problem (2 for 2D, 3 for 3D).
      %
      % Outputs:
      %   This function does not return any outputs but modifies properties of
      %   the `obj` instance.
      %
      % Uses:
      %   obj.nmixcell - (int32) The number of mixed cells in the problem. Used to
      %                  determine whether to allocate or clear storage.
      %
      % Modifies:
      %   obj.mix_mpoly_to_mix
      %      (int array) Mapping from mixed cell index to the mixed-cell-with-
      %                  multipolygons index.
      %
      %   obj.nnode_in_mixcell
      %      (int array) The number of nodes in each mixed cell.
      %
      %   obj.coords_in_mixcell
      %      (double array) Coordinates of nodes in each mixed
      %                  cell, with each cell containing a double array.
      %
      %   obj.nnode_for_minterface_in_mixcell
      %      (cell array) The number of nodes for each material interface in mixed
      %                  cells, with each cell containing an int array.
      %
      %   obj.nodes_for_minterface_in_mixcell
      %      (cell array) The nodes for each material interface in mixed cells,
      %                  with each cell containing an int array.
      %
      %   obj.nnode_for_mpoly_in_mixcell
      %      (2D cell array) The number of nodes for each material
      %                  polygon in mixed cells.
      %
      %   obj.nodes_for_mpoly_in_mixcell
      %      (2D cell array) The nodes for each material polygon in mixed cells.
      %
      %   obj.nface_for_mpoly_in_mixcell
      %      (2D cell array) The number of faces for each material polyhedron
      %                  in mixed cells.
      %
      %   obj.nnode_for_face_ea_mpoly_in_mixcell
      %      (3D cell array) The number of nodes for each face of a material
      %                  polyhedron in mixed cells.
      %
      %   obj.nodelist_for_face_ea_mpoly_in_mixcell
      %      (3D cell array) The nodes for each face of a material polyhedron
      %                  in each mixed cell.
      %
      % Original C declaration:
      % void set_mpoly_storage(int dim)

      % Clean up existing multipolygon storage
      obj.clean_mpoly();

      if obj.nmixcell > 0
        obj.mix_mpoly_to_mix = zeros(obj.nmixcell, 1);
        obj.nnode_in_mixcell = zeros(obj.nmixcell, 1);
        obj.coords_in_mixcell = cell(obj.nmixcell, 1);
        obj.nnode_for_minterface_in_mixcell = cell(obj.nmixcell, 1);
        obj.nodes_for_minterface_in_mixcell = cell(obj.nmixcell, 1);

        if dim == 2
          obj.nnode_for_mpoly_in_mixcell = cell(obj.nmixcell, 1);
          obj.nodes_for_mpoly_in_mixcell = cell(obj.nmixcell, 1);

          for mx = 1:obj.nmixcell
            obj.coords_in_mixcell{mx} = [];
            obj.nnode_for_minterface_in_mixcell{mx} = [];
            obj.nodes_for_minterface_in_mixcell{mx} = [];
            obj.nnode_for_mpoly_in_mixcell{mx} = [];
            obj.nodes_for_mpoly_in_mixcell{mx} = {};
          end
        elseif dim == 3
          obj.nface_for_mpoly_in_mixcell = cell(obj.nmixcell, 1);
          obj.nnode_for_face_ea_mpoly_in_mixcell = cell(obj.nmixcell, 1);
          obj.nodelist_for_face_ea_mpoly_in_mixcell = cell(obj.nmixcell, 1);

          for mx = 1:obj.nmixcell
            obj.coords_in_mixcell{mx} = [];
            obj.nnode_for_minterface_in_mixcell{mx} = [];
            obj.nodes_for_minterface_in_mixcell{mx} = [];
            obj.nface_for_mpoly_in_mixcell{mx} = {};
            obj.nnode_for_face_ea_mpoly_in_mixcell{mx} = {};
            obj.nodelist_for_face_ea_mpoly_in_mixcell{mx} = {};
          end
        end
      else
        obj.mix_mpoly_to_mix = [];
        obj.nnode_in_mixcell = [];
        obj.coords_in_mixcell = {};
        obj.nnode_for_minterface_in_mixcell = {};
        obj.nodes_for_minterface_in_mixcell = {};
        obj.nnode_for_mpoly_in_mixcell = [];
        obj.nodes_for_mpoly_in_mixcell = {};
        obj.nface_for_mpoly_in_mixcell = {};
        obj.nnode_for_face_ea_mpoly_in_mixcell = {};
        obj.nodelist_for_face_ea_mpoly_in_mixcell = {};
      end
    end

    function set_2dmesh_mat(obj, mesh, nreg, is_solid, gamma_ea_mat, ...
        matids_ea_reg, reg_shape, rho_ea_reg, pres_ea_reg, ...
        ei_ea_reg, v_ea_reg)
      % Function: set_2dmesh_mat
      %
      % This function sets up the material properties for every cell on a
      % 2D mesh.  It takes, as an input, a number of regions described by
      % their respective geometric shapes, and material assignments for
      % each of the regions: material ID, density, pressure, internal
      % energy, velocity etc. It then proceeds with calculating the volume
      % fractions of each material within every cell.
      %
      % Inputs:
      % - obj (c_mat): The current c_mat object instance.
      % - mesh (c_mesh): An instance of the c_mesh class.
      % - nreg (int): Number of regions.
      % - is_solid (int array): Indicator array for solid materials.
      % - matids_ea_reg (int array): Mapping regions to material IDs.
      % - gamma_ea_mat (double array): Gamma values for each material.
      % - reg_shape (region_shape array): Geometric shapes defining regions.
      % - rho_ea_reg (double array): Density for each region.
      % - pres_ea_reg (double array): Pressure for each region.
      % - ei_ea_reg (double array): Internal energy density for each region.
      % - v_ea_reg (cell): Cell array containing velocity arrays for each
      %   region.
      %
      % Outputs:
      % - None (modifies variables within obj and mesh).
      %
      % Uses:
      % - mesh.ncell_prob (int array): Number of cells in each dimension.
      % - mesh.dx (double array): Cell sizes in each dimension.
      % - mesh.xl_prob (double array): Lower bounds of the simulation domain.
      % - mesh.nbdry_prob (int): Number of boundary layers.
      % - mesh.nmat_mesh (int): Number of materials in the mesh.
      % - mesh.matids_mesh (int array): Material IDs in the mesh.
      %
      % Modifies:
      % - obj.nmixcell_int (int): The number of mixed cells in interior cells.
      % - obj.nmixcell (int): The number of mixed cells in the problem.
      % - obj.nmat_in_mixcell (int array): Number of materials in each mixed
      %   cell.
      % - obj.ijk_in_mixcell (int array): The cell indices of mixed cells.
      % - obj.matids_in_mixcell (cell array of int arrays): Material IDs
      %   for each material in each mixed cell.
      % - obj.vf_in_mixcell (cell array of double arrays): Volume fractions
      %   of each material in each mixed cell.
      % - obj.rho_in_mixcell (cell array of double arrays): Density for each
      %   mixed cell.
      % - obj.pres_in_mixcell (cell array of double arrays): Pressure for each
      %   mixed cell.
      % - obj.ei_in_mixcell (cell array of double arrays): Internal energy
      %   density for each mixed cell.
      % - mesh.nmat_for_2dcell (int array): Number of materials in each 2D cell.
      % - mesh.matid_for_2dcell (int array): Material IDs for each 2D cell.
      % - mesh.mixcell_for_2dcell (int array): Mixed cell indices for 2D cells.
      % - mesh.rho_for_2dcell (double array): Density for each 2D cell.
      % - mesh.ei_for_2dcell (double array): Internal energy density for each
      %   2D cell.
      % - mesh.pres_for_2dcell (double array): Pressure for each 2D cell.
      % - vel_for_2dcell (double array): Velocity for each 2D cell.
      %
      % Original C declaration:
      % void set_2dmesh_mat(int *ncell, double *xl_prob, double *dx,
      %                     int nbdry, int nmat, int *matids_ea_mat,
      %                     int *is_solid, double *gamma_ea_mat, int nreg,
      %                     int *matids_ea_reg, Region_Shape *reg_shape,
      %                     double *rho_ea_reg, double *pres_ea_reg,
      %                     double *ei_ea_reg, double **v_ea_reg,
      %                     int **nmat_for_2dcell, int **matid_for_2dcell,
      %                     int **mixcell_for_2dcell, double **rho_for_2dcell,
      %                     double **ei_for_2dcell, double **pres_for_2dcell,
      %                     double ***vel_for_2dcell)

      % Initialize variables
      nmix = 0;
      nmix_mx = 1.0;  % Estimate the number of mixed cells
      dim = 2;
      ncell_ext = zeros(1, dim);
      dx = mesh.dx(1:dim);

      for i = 1:dim
        ncell_ext(i) = mesh.ncell_prob(i) + mesh.nbdry_prob * 2;
        nmix_mx = nmix_mx * (0.1 * ncell_ext(i) + 1);
      end
      nmix_mx = round(nmix_mx);

      % Estimate space for mixed cells
      lsize_mx = nmix_mx * nreg + 1;
      vfs_list = zeros(lsize_mx, 1);
      rho_list = zeros(lsize_mx, 1);
      pres_list = zeros(lsize_mx, 1);
      ei_list = zeros(lsize_mx, 1);
      ids_list = zeros(lsize_mx, 1);

      lsize = 0; % The actual length in vfs_list and ids_list

      % Estimate the number of mix cells
      tmp_ijk_in_mixcell = zeros(nmix_mx + 1, dim);
      obj.nmat_in_mixcell = zeros(nmix_mx + 1, 1);

      xl(2) = mesh.xl_prob(2) - mesh.dx(2) * mesh.nbdry_prob;
      for j = 1:ncell_ext(2)

        xl(1) = mesh.xl_prob(1) - mesh.dx(1) * mesh.nbdry_prob;
        for i = 1:ncell_ext(1)
          [nmat_cell, matids_cell, rho, pres, ei, vel, ...
            vf_ea_mat, rho_ea_mat, pres_ea_mat, ei_ea_mat] ...
            = obj.set_cell_mat(dim, xl, dx, nreg, is_solid, ...
            matids_ea_reg, gamma_ea_mat, reg_shape, ...
            rho_ea_reg, pres_ea_reg, ei_ea_reg, ...
            v_ea_reg);
          nm = 0;
          vfsum = 0.0;

          for m = 1:nmat_cell
            if vf_ea_mat(m) < obj.vfmin
              vf_ea_mat(m) = 0.0;
            else
              vfsum = vfsum + vf_ea_mat(m);
              nm = nm + 1;
            end
          end

          vfsuminv = 1.0 / vfsum;
          vf_ea_mat(1:nmat_cell) = vf_ea_mat(1:nmat_cell) * vfsuminv;

          mesh.nmat_for_2dcell(j, i) = nm;
          mesh.matid_for_2dcell(j, i) = matids_cell(1);

          if nm > 1
            mesh.mixcell_for_2dcell(j, i) = nmix + 1;
            % If the number of mixed cells exceeded the
            % allocated limit:
            if nmix >= nmix_mx
              nmix_mx = nmix_mx + floor(0.1 * nmix_mx + 10);
              obj.nmat_in_mixcell = [obj.nmat_in_mixcell; zeros(nmix_mx - length(obj.nmat_in_mixcell), 1)];
              tmp_ijk_in_mixcell = [tmp_ijk_in_mixcell; zeros(nmix_mx - size(tmp_ijk_in_mixcell, 1), dim)];
            end

            obj.nmat_in_mixcell(nmix + 1) = nm;
            tmp_ijk_in_mixcell(nmix + 1, :) = [i, j];

            if lsize + nm >= lsize_mx
              lsize_mx = lsize_mx + floor(0.1 * lsize_mx + 10 * nm);
              vfs_list = [vfs_list; zeros(lsize_mx - length(vfs_list), 1)];
              rho_list = [rho_list; zeros(lsize_mx - length(rho_list), 1)];
              pres_list = [pres_list; zeros(lsize_mx - length(pres_list), 1)];
              ei_list = [ei_list; zeros(lsize_mx - length(ei_list), 1)];
              ids_list = [ids_list; zeros(lsize_mx - length(ids_list), 1)];
            end

            for m = 1:nmat_cell
              if vf_ea_mat(m) > 0.0
                lsize = lsize + 1;
                vfs_list(lsize) = vf_ea_mat(m);
                rho_list(lsize) = rho_ea_mat(m);
                pres_list(lsize) = pres_ea_mat(m);
                ei_list(lsize) = ei_ea_mat(m);
                ids_list(lsize) = matids_cell(m);
              end
            end
            nmix = nmix + 1;
          end

          mesh.nmat_for_2dcell(j, i) = nmat_cell;
          mesh.matid_for_2dcell(j, i) = matids_cell(1);
          mesh.rho_for_2dcell(j, i) = rho;
          mesh.ei_for_2dcell(j, i) = ei;
          mesh.pres_for_2dcell(j, i) = pres;
          mesh.vel_for_2dcell(j, i, :) = vel;

          xl(1) = xl(1) + mesh.dx(1);
        end
        xl(2) = xl(2) + mesh.dx(2);
      end

      mesh.mesh_sound_speed(mesh.nmat_mesh, mesh.matids_mesh, is_solid, gamma_ea_mat)

      if nmix == 0
        % If there are no mixed cells, set the fields to empty arrays
        obj.ijk_in_mixcell = [];
        obj.nmat_in_mixcell = [];
      else
        % Resize nmat_in_mixcell to match the actual number of mixed cells (nmix)
        obj.nmat_in_mixcell = obj.nmat_in_mixcell(1:nmix);

        % Resize ijk_in_mixcell to match the actual number of mixed cells (nmix)
        obj.ijk_in_mixcell = zeros(nmix, dim);
        obj.ijk_in_mixcell(1:nmix,1:dim) = tmp_ijk_in_mixcell(1:nmix,1:dim);

        % Allocate memory for the mixed cell data using cell arrays
        obj.vf_in_mixcell = cell(nmix, 1);
        obj.rho_in_mixcell = cell(nmix, 1);
        obj.pres_in_mixcell = cell(nmix, 1);
        obj.ei_in_mixcell = cell(nmix, 1);
        obj.matids_in_mixcell = cell(nmix, 1);

        % Assign variable-sized arrays to each cell
        start_idx = 1;
        for s = 1:nmix
          end_idx = start_idx + obj.nmat_in_mixcell(s) - 1;

          % Populate the vf_in_mixcell, rho_in_mixcell, pres_in_mixcell, ei_in_mixcell, matids_in_mixcell arrays
          obj.vf_in_mixcell{s} = vfs_list(start_idx:end_idx);
          obj.rho_in_mixcell{s} = rho_list(start_idx:end_idx);
          obj.pres_in_mixcell{s} = pres_list(start_idx:end_idx);
          obj.ei_in_mixcell{s} = ei_list(start_idx:end_idx);
          obj.matids_in_mixcell{s} = ids_list(start_idx:end_idx);

          start_idx = end_idx + 1;
        end
      end

      obj.nmixcell_int = nmix;
      obj.nmixcell = nmix;

    end

    function set_3dmesh_mat(obj, mesh, nreg, is_solid, gamma_ea_mat, ...
        matids_ea_reg, reg_shape, rho_ea_reg, pres_ea_reg, ...
        ei_ea_reg, v_ea_reg)
      % Function: set_3dmesh_mat
      %
      % This function sets up the material properties for every cell on a
      % 3D mesh.  It takes, as an input, a number of regions described by
      % their respective geometric shapes, and material assignments for
      % each of the regions: material ID, density, pressure, internal
      % energy, velocity etc. It then proceeds with calculating the volume
      % fractions of each material within every cell.
      %
      % Inputs:
      % - obj (c_mat): The current c_mat object instance.
      % - mesh (c_mesh): An instance of the c_mesh class.
      % - nreg (int): Number of regions.
      % - is_solid (int array): Indicator array for solid materials.
      % - matids_ea_reg (int array): Mapping regions to material IDs.
      % - gamma_ea_mat (double array): Gamma values for each material.
      % - reg_shape (region_shape array): Geometric shapes defining regions.
      % - rho_ea_reg (double array): Density for each region.
      % - pres_ea_reg (double array): Pressure for each region.
      % - ei_ea_reg (double array): Internal energy density for each region.
      % - v_ea_reg (cell): Cell array containing velocity arrays for each
      %   region.
      %
      % Outputs:
      % - None (modifies variables within obj and mesh).
      %
      % Uses:
      % - mesh.ncell_prob (int array): Number of cells in each dimension.
      % - mesh.dx (double array): Cell sizes in each dimension.
      % - mesh.xl_prob (double array): Lower bounds of the simulation domain.
      % - mesh.nbdry_prob (int): Number of boundary layers.
      % - mesh.nmat_mesh (int): Number of materials in the mesh.
      % - mesh.matids_mesh (int array): Material IDs in the mesh.
      %
      % Modifies:
      % - obj.nmixcell_int (int): The number of mixed cells in interior cells.
      % - obj.nmixcell (int): The number of mixed cells in the problem.
      % - obj.nmat_in_mixcell (int array): Number of materials in each mixed
      %   cell.
      % - obj.ijk_in_mixcell (int array): The cell indices of mixed cells.
      % - obj.matids_in_mixcell (cell array of int arrays): Material IDs
      %   for each material in each mixed cell.
      % - obj.vf_in_mixcell (cell array of double arrays): Volume fractions
      %   of each material in each mixed cell.
      % - obj.rho_in_mixcell (cell array of double arrays): Density for each
      %   mixed cell.
      % - obj.pres_in_mixcell (cell array of double arrays): Pressure for each
      %   mixed cell.
      % - obj.ei_in_mixcell (cell array of double arrays): Internal energy
      %   density for each mixed cell.
      % - mesh.nmat_for_3dcell (int array): Number of materials in each 3D cell.
      % - mesh.matid_for_3dcell (int array): Material IDs for each 3D cell.
      % - mesh.mixcell_for_3dcell (int array): Mixed cell indices for 3D cells.
      % - mesh.rho_for_3dcell (double array): Density for each 3D cell.
      % - mesh.ei_for_3dcell (double array): Internal energy density for each
      %   3D cell.
      % - mesh.pres_for_3dcell (double array): Pressure for each 3D cell.
      % - vel_for_3dcell (double array): Velocity for each 3D cell.
      %
      % Original C declaration:
      % void set_3dmesh_mat(int *ncell, double *xl_prob, double *dx,
      %                     int nbdry, int nmat, int *matids_ea_mat,
      %                     int *is_solid, double *gamma_ea_mat, int nreg,
      %                     int *matids, Region_Shape *reg_shape,
      %                     double *rho_ea_reg, double *pres_ea_reg,
      %                     double *ei_ea_reg, double **v_ea_reg,
      %                     int ***nmat_for_3dcell, int ***matid_for_3dcell,
      %                     int ***mixcell_for_3dcell, double ***rho_for_3dcell,
      %                     double ***ei_for_3dcell, double ***pres_for_3dcell,
      %                     double ****vel_for_3dcell);

      % Initialize variables
      nmix = 0;
      nmix_mx = 1.0;  % Estimate the number of mixed cells
      dim = 3;
      ncell_ext = zeros(1, dim);
      dx = mesh.dx(1:dim);

      for i = 1:dim
        ncell_ext(i) = mesh.ncell_prob(i) + mesh.nbdry_prob * 2;
        nmix_mx = nmix_mx * (0.1 * ncell_ext(i) + 1);
      end
      nmix_mx = round(nmix_mx);

      % Estimate space for mixed cells
      lsize_mx = nmix_mx * nreg + 1;
      vfs_list = zeros(lsize_mx, 1);
      rho_list = zeros(lsize_mx, 1);
      pres_list = zeros(lsize_mx, 1);
      ei_list = zeros(lsize_mx, 1);
      ids_list = zeros(lsize_mx, 1);

      lsize = 0; % The actual length in vfs_list and ids_list

      % Estimate the number of mix cells
      tmp_ijk_in_mixcell = zeros(nmix_mx + 1, dim);
      obj.nmat_in_mixcell = zeros(nmix_mx + 1, 1);

      xl(3) = mesh.xl_prob(3) - mesh.dx(3) * mesh.nbdry_prob;
      for k = 1:ncell_ext(3)

        xl(2) = mesh.xl_prob(2) - mesh.dx(2) * mesh.nbdry_prob;
        for j = 1:ncell_ext(2)

          xl(1) = mesh.xl_prob(1) - mesh.dx(1) * mesh.nbdry_prob;
          for i = 1:ncell_ext(1)
            [nmat_cell, matids_cell, rho, pres, ei, vel, ...
              vf_ea_mat, rho_ea_mat, pres_ea_mat, ei_ea_mat] ...
              = obj.set_cell_mat(dim, xl, dx, nreg, is_solid, ...
              matids_ea_reg, gamma_ea_mat, reg_shape, ...
              rho_ea_reg, pres_ea_reg, ei_ea_reg, ...
              v_ea_reg);
            nm = 0;
            vfsum = 0.0;

            for m = 1:nmat_cell
              if vf_ea_mat(m) < obj.vfmin
                vf_ea_mat(m) = 0.0;
              else
                vfsum = vfsum + vf_ea_mat(m);
                nm = nm + 1;
              end
            end

            vfsuminv = 1.0 / vfsum;
            vf_ea_mat(1:nmat_cell) = vf_ea_mat(1:nmat_cell) * vfsuminv;

            mesh.nmat_for_3dcell(j, i) = nm;
            for m = 1:nmat_cell
              if vf_ea_mat(m) > 0.0
                mesh.matid_for_3dcell(k, j, i) = matids_cell(m);
                break;
              end
            end

            if nm > 1
              mesh.mixcell_for_2dcell(j, i) = nmix + 1;
              % If the number of mixed cells exceeded the
              % allocated limit:
              if nmix >= nmix_mx
                nmix_mx = nmix_mx + floor(0.1 * nmix_mx + 10);
                obj.nmat_in_mixcell = [obj.nmat_in_mixcell; zeros(nmix_mx - length(obj.nmat_in_mixcell), 1)];
                tmp_ijk_in_mixcell = [tmp_ijk_in_mixcell; zeros(nmix_mx - size(tmp_ijk_in_mixcell, 1), dim)];
              end

              obj.nmat_in_mixcell(nmix + 1) = nm;
              tmp_ijk_in_mixcell(nmix + 1, :) = [i, j, k];

              if lsize + nm >= lsize_mx
                lsize_mx = lsize_mx + floor(0.1 * lsize_mx + 10 * nm);
                vfs_list = [vfs_list; zeros(lsize_mx - length(vfs_list), 1)];
                rho_list = [rho_list; zeros(lsize_mx - length(rho_list), 1)];
                pres_list = [pres_list; zeros(lsize_mx - length(pres_list), 1)];
                ei_list = [ei_list; zeros(lsize_mx - length(ei_list), 1)];
                ids_list = [ids_list; zeros(lsize_mx - length(ids_list), 1)];
              end

              for m = 1:nmat_cell
                if vf_ea_mat(m) > 0.0
                  lsize = lsize + 1;
                  vfs_list(lsize) = vf_ea_mat(m);
                  rho_list(lsize) = rho_ea_mat(m);
                  pres_list(lsize) = pres_ea_mat(m);
                  ei_list(lsize) = ei_ea_mat(m);
                  ids_list(lsize) = matids_cell(m);
                end
              end
              nmix = nmix + 1;
            end

            mesh.nmat_for_3dcell(k, j, i) = nmat_cell;
            mesh.matid_for_3dcell(k, j, i) = matids_cell(1);
            mesh.rho_for_3dcell(k, j, i) = rho;
            mesh.ei_for_3dcell(k, j, i) = ei;
            mesh.pres_for_3dcell(k, j, i) = pres;
            mesh.vel_for_3dcell(k, j, i, :) = vel;

            xl(1) = xl(1) + mesh.dx(1);
          end
          xl(2) = xl(2) + mesh.dx(2);
        end
        xl(3) = xl(3) + mesh.dx(3);
      end

      mesh.mesh_sound_speed(mesh.nmat_mesh, mesh.matids_mesh, is_solid, gamma_ea_mat)

      if nmix == 0
        % If there are no mixed cells, set the fields to empty arrays
        obj.ijk_in_mixcell = [];
        obj.nmat_in_mixcell = [];
      else
        % Resize nmat_in_mixcell to match the actual number of mixed cells (nmix)
        obj.nmat_in_mixcell = obj.nmat_in_mixcell(1:nmix);

        % Resize ijk_in_mixcell to match the actual number of mixed cells (nmix)
        obj.ijk_in_mixcell = zeros(nmix, dim);
        obj.ijk_in_mixcell(1:nmix,1:dim) = tmp_ijk_in_mixcell(1:nmix,1:dim);

        % Allocate memory for the mixed cell data using cell arrays
        obj.vf_in_mixcell = cell(nmix, 1);
        obj.rho_in_mixcell = cell(nmix, 1);
        obj.pres_in_mixcell = cell(nmix, 1);
        obj.ei_in_mixcell = cell(nmix, 1);
        obj.matids_in_mixcell = cell(nmix, 1);

        % Assign variable-sized arrays to each cell
        start_idx = 1;
        for s = 1:nmix
          end_idx = start_idx + obj.nmat_in_mixcell(s) - 1;

          % Populate the vf_in_mixcell, rho_in_mixcell, pres_in_mixcell, ei_in_mixcell, matids_in_mixcell arrays
          obj.vf_in_mixcell{s} = vfs_list(start_idx:end_idx);
          obj.rho_in_mixcell{s} = rho_list(start_idx:end_idx);
          obj.pres_in_mixcell{s} = pres_list(start_idx:end_idx);
          obj.ei_in_mixcell{s} = ei_list(start_idx:end_idx);
          obj.matids_in_mixcell{s} = ids_list(start_idx:end_idx);

          start_idx = end_idx + 1;
        end
      end

      obj.nmixcell_int = nmix;
      obj.nmixcell = nmix;

    end

    function [nmat_cell, matids_cell, rho, pres, ei, vel, ...
        vf_ea_mat, rho_ea_mat, pres_ea_mat, ei_ea_mat] = ...
        set_cell_mat(obj, dim, xl, dx, nreg, is_solid, ...
        matids_ea_reg, gamma_ea_mat, reg_shape, ...
        rho_ea_reg, pres_ea_reg, ei_ea_reg, ...
        v_ea_reg)
      % Function: set_cell_mat
      %
      % Calculates the material properties for a given cell in the mesh,
      % given the regions and materials defined in the input.  Cell
      % properties include the total number of materials in the cell,
      % their volume fractions, densities, pressures, internal energies,
      % and velocities.
      %
      % Inputs:
      % - obj (c_mat): The current c_mat object instance.
      % - dim (double): Dimensionality of the problem.
      % - xl (double array): Lower coordinates of the cell [size: dim, 1].
      % - dx (double array): Cell dimensions [size: dim, 1].
      % - nreg (int): Number of regions.
      % - is_solid (int array): Indicator array for solid materials [size: nmat].
      % - matids_ea_reg (int array): Mapping regions to material IDs.
      % - gamma_ea_mat (double array): Gamma values for each material [size: nmat].
      % - reg_shape (region_shape array): Array of region_shape objects defining the
      %   geometry of each region.
      % - rho_ea_reg (double array): Densities for each region [size: nreg, 1].
      % - pres_ea_reg (double array): Pressures for each region [size: nreg, 1].
      % - ei_ea_reg (double array): Internal energy densities for each region [size:
      %   nreg, 1].
      % - v_ea_reg (cell): Cell array of velocity arrays for each region [size: nreg,
      %   1], each cell containing a double array [size: dim, 1].
      %
      % Outputs:
      % - nmat_cell (double): Number of materials in the cell.
      % - matids_cell (double array): Material IDs in the cell [size: nmat_cell, 1].
      % - rho (double): Total density in the cell.
      % - pres (double): Total pressure in the cell.
      % - ei (double): Total internal energy in the cell.
      % - vel (double array): Velocity in the cell [size: dim, 1].
      % - vf_ea_mat (double array): Volume fractions for each material in the cell
      %   [size: nmat_cell, 1].
      % - rho_ea_mat (double array): Densities for each material in the cell [size:
      %   nmat_cell, 1].
      % - pres_ea_mat (double array): Pressures for each material in the cell [size:
      %   nmat_cell, 1].
      % - ei_ea_mat (double array): Internal energies for each material in the cell
      %   [size: nmat_cell, 1].
      %
      % Uses:
      % - obj.vfmin (double): Smallest volume fraction, otherwise considered zero.
      %
      % Original C declaration:
      % void set_cell_mat(int dim, double *xl, double *dx,
      %                   int nmat, int *matids_ea_mat, int *is_solid,
      %                   double *gamma_ea_mat, int nreg, int *matids_ea_reg,
      %                   Region_Shape *reg_shape, double *rho_ea_reg,
      %                   double *pres_ea_reg, double *ei_ea_reg,
      %                   double **v_ea_reg, int *nmat_cell, int *matids_cell,
      %                   double *vf_ea_mat, double *rho_ea_mat,
      %                   double *pres_ea_mat, double *ei_ea_mat,
      %                   double *rho_cell, double *pres_cell,
      %                   double *ei_cell, double *v_cell)

      % Constants and allocations
      global tiny;
      global geom eos;

      geop = 1;       % Cartesian coordinates
      ifinquiry = 0;  % to get vol fraction, not just mixed or clean
      mixed = 0;

      nmat = obj.nmat_prob;
      vf_ea_mat = zeros(nmat, 1);
      rho_ea_mat = zeros(nmat, 1);
      pres_ea_mat = zeros(nmat, 1);
      ei_ea_mat = zeros(nmat, 1);
      matids_cell = zeros(nmat, 1);

      % Initialize variables
      cell_vol = prod(dx);
      vf_ea_mat(1) = 1.0;

      for r = 2:nreg
        vf_ea_mat(r) = 0.0;
        switch reg_shape(r).type
          case shape_type.shape_sphere
            rad1 = reg_shape(r).parameters(1);
            ctr1 = reg_shape(r).parameters(2:end);
            [~, vol] = geom.gsph_rec(ifinquiry, geop, dim, ctr1, rad1, xl, dx);

          case shape_type.shape_quad
            nn = 4;
            [xmin, xmax] = obj.find_min_max(dim, nn, reg_shape(r).parameters);
            [~, vol] = geom.poly2d_rec(ifinquiry, geop, cell_vol, dim, nn, reg_shape(r).parameters, xmin, xmax, xl, dx);

          case shape_type.shape_rectangular
            xl_reg = reg_shape(r).parameters(1:dim);
            xr_reg = reg_shape(r).parameters(dim+1:end);
            [~, vol] = geom.rec_rec(ifinquiry, geop, cell_vol, dim, xl_reg, xr_reg, xl, dx);

          case shape_type.shape_cylinder
            rad1 = reg_shape(r).parameters(1);
            ctr1 = reg_shape(r).parameters(2:dim+1);
            rad2 = reg_shape(r).parameters(dim+2);
            ctr2 = reg_shape(r).parameters(dim+3:end);
            ifcyl = find(abs(ctr1 - ctr2) > 0.0001 * rad1, 1);
            if isempty(ifcyl)
              ifcyl = 0;
            else
              ifcyl = ifcyl + 1;
            end
            [~, vol] = geom.gconic_rec(ifinquiry, geop, ifcyl, dim, rad1, rad2, ctr1, ctr2, xl, dx);

          otherwise
            error('ERROR: reg type not yet implemented');
        end
        vf = vol / cell_vol;
        vf = max(min(1.0, vf), 0.0);
        vf_left = vf;
        vf_ea_mat(r) = vf;
        previous_reg = r - 1;
        while vf_left > 0.0 && previous_reg >= 1
          if vf_left <= vf_ea_mat(previous_reg)
            vf_ea_mat(previous_reg) = vf_ea_mat(previous_reg) - vf_left;
            vf_left = 0.0;
          else
            vf_left = vf_left - vf_ea_mat(previous_reg);
            vf_ea_mat(previous_reg) = 0.0;
          end
          previous_reg = previous_reg - 1;
        end
      end

      % Calculate properties of the cell
      rho = 0.0;
      ei = 0.0;
      pres = 0.0;
      mass = 0.0;
      vel = zeros(dim, 1);

      for r = 1:nreg
        rho = rho + vf_ea_mat(r) * rho_ea_reg(r);
        ei = ei + vf_ea_mat(r) * ei_ea_reg(r);
        pres = pres + vf_ea_mat(r) * pres_ea_reg(r);
        mass = mass + vf_ea_mat(r) * rho_ea_reg(r);
        for i = 1:dim
          vel(i) = vel(i) + vf_ea_mat(r)*rho_ea_reg(r)*v_ea_reg{r}(i);
        end
      end

      % Recover velocity
      massinv = 1.0 / (mass + tiny);
      vel = vel * massinv;

      % Consolidate materials
      nmat_cell = 0;
      for r = 1:nreg
        mass = 0.0;
        ener = 0.0;
        if vf_ea_mat(r) > 0.0
          matid = matids_ea_reg(r);
          vf = vf_ea_mat(r);
          mass = mass + vf_ea_mat(r) * rho_ea_reg(r);
          ener = ener + vf_ea_mat(r) * ei_ea_reg(r);

          for i = r+1:nreg
            if matids_ea_reg(i) == matid
              vf = vf + vf_ea_mat(i);
              mass = mass + vf_ea_mat(i) * rho_ea_reg(i);
              ener = ener + vf_ea_mat(i) * ei_ea_reg(i);
              vf_ea_mat(i) = 0.0;
            end
          end
          rho_ea_mat(nmat_cell + 1) = mass / vf;
          ei_ea_mat(nmat_cell + 1) = ener / vf;
          if ~is_solid(matid)
            pres_ea_mat(nmat_cell + 1) ...
              = (gamma_ea_mat(matid) - 1.0) ...
              * ei_ea_mat(nmat_cell + 1);
          else
            pres_ea_mat(nmat_cell + 1) ...
              = eos.p_mie_gruneisen(rho_ea_mat(nmat_cell + 1),...
              ei_ea_mat(nmat_cell + 1));
          end
          matids_cell(nmat_cell + 1) = matid;
          vf_ea_mat(nmat_cell + 1) = vf;
          nmat_cell = nmat_cell + 1;
        end
      end
    end

    function cal_mixcell_zgrad2d(obj, mesh)
      % cal_mixcell_zgrad2d - Calculate the gradient of volume fractions for each
      %                       material in mixed 2D cells using properties from
      %                       c_mat and c_mesh objects.
      %
      % This function calculates the gradient of volume fractions for each material
      % in mixed cells within a 2D grid. The gradients are computed for each material
      % in each mixed zone based on the surrounding cells. The function leverages
      % properties from both the `c_mat` and `c_mesh` classes and calls
      % `cal_cell_zgrad2d` to compute the gradient for each individual cell.
      %
      % Inputs:
      %   obj                   - (c_mat object) The instance of the c_mat class
      %                           containing material data and properties such as
      %                           material IDs, volume fractions, and mixed cell data.
      %   mesh                  - (c_mesh object) The instance of the c_mesh class
      %                           containing mesh-related properties such as cell
      %                           counts, boundary layers, and material information.
      %
      % Outputs:
      %   This function does not return any outputs but modifies the `vfgrad_mixcell`
      %   property of the `c_mat` class.
      %
      % Uses:
      %   obj.nmixcell             - (int32) The number of mixed cells in the problem.
      %   obj.nmat_in_mixcell      - (int array) The number of materials in each mixed zone.
      %   obj.matids_in_mixcell    - (cell array of int arrays) The material IDs in each
      %                              mixed zone.
      %   obj.vf_in_mixcell        - (cell array of double arrays) The volume fractions
      %                              of each material in each mixed zone.
      %   obj.ijk_in_mixcell       - (array of int arrays) The indices of each mixed
      %                              zone in the grid.
      %   mesh.ncell_prob          - (int array) The number of cells in the x and y
      %                              directions, size: [2, 1].
      %   mesh.nbdry_prob          - (int) The number of boundary cells.
      %   mesh.xl_prob             - (double array) The lower bounds of the problem
      %                              domain in the x and y directions, size: [2, 1].
      %   mesh.xr_prob             - (double array) The upper bounds of the problem
      %                              domain in the x and y directions, size: [2, 1].
      %   mesh.nmat_for_2dcell     - (cell array of int arrays) The number of materials
      %                              in each 2D cell.
      %   mesh.matid_for_2dcell    - (cell array of int arrays) The material IDs in each
      %                              2D cell.
      %   mesh.mixcell_for_2dcell  - (cell array of int arrays) The mixed cell indices
      %                              for each 2D cell.
      %
      % Modifies:
      %   obj.vfgrad_mixcell       - (cell array of double arrays) The gradient of the
      %                              volume fractions for each material in each mixed zone.
      %
      % Original C declaration:
      % void cal_mixcell_zgrad2d(int *ncell, int nbdry, double *xl_prob, double *dx,
      %                          int **nmat_for_2dcell, int **matid_for_2dcell, int **mixcell_for_2dcell,
      %                          int nmixzone, int *nmat_for_mzone, int **matid_for_mzone, double **vf_for_mzone,
      %                          int **ijk_mzone, double ***vfgrad_for_mzone)

      % Calculate the last cell index accounting for boundary layers
      dim = 2;
      ncell_last = mesh.ncell_prob + 2 * mesh.nbdry_prob;

      % Loop over all mixed zones
      for mix = 1:obj.nmixcell
        nm = obj.nmat_in_mixcell(mix);
        matids = obj.matids_in_mixcell{mix};
        vfs = obj.vf_in_mixcell{mix};
        ijk = obj.ijk_in_mixcell(mix,1:dim);

        for m = 1:nm
          i0 = ijk(1);
          j0 = ijk(2);

          % Skip boundary cells
          if i0 == 1 || i0 == ncell_last(1) ...
              || j0 == 1 || j0 == ncell_last(2)
            continue;
          end

          % Initialize var and calculate the gradient
          matid = matids(m);
          var = zeros(3, 3);
          var(2, 2) = vfs(m);

          for j = -1:1
            jcell = j0 + j;
            j3 = j + 2;  % index in var[3][3]
            for i = -1:1
              if i == 0 && j == 0
                continue;
              end
              icell = i0 + i;
              i3 = i + 2;  % index in var[3][3]

              var(j3, i3) = 0.0;
              nm_cell = mesh.nmat_for_2dcell(jcell, icell);
              matid_cell = mesh.matid_for_2dcell(jcell, icell);
              if nm_cell == 1
                if matid_cell == matid
                  var(j3, i3) = 1.0;
                end
              else
                mix_nb = mesh.mixcell_for_2dcell(jcell, icell);
                nm_nb = obj.nmat_in_mixcell(mix_nb);
                assert(nm_nb == nm_cell);
                matids_nb = obj.matids_in_mixcell{mix_nb};
                vfs_nb = obj.vf_in_mixcell{mix_nb};
                for m_nb = 1:nm_nb
                  if matid == matids_nb(m_nb)
                    var(j3, i3) = vfs_nb(m_nb);
                    break;
                  end
                end
              end
            end
          end

          % Calculate the gradient for the current cell
          gradv = mesh.cal_cell_zgrad2d(mesh.xr_prob - mesh.xl_prob, var);
          for i = 1:dim
            obj.vfgrad_mixcell{mix}(m, i) = gradv(i);
          end
        end
      end
    end

    function cal_mixcell_zgrad3d(obj, mesh)
      % cal_mixcell_zgrad3d - Calculate the gradient of volume fractions for each
      %                       material in mixed 3D cells using properties from
      %                       c_mat and c_mesh objects.
      %
      % This function calculates the gradient of volume fractions for each material
      % in mixed cells within a 3D grid. The gradients are computed for each material
      % in each mixed zone based on the surrounding cells. The function leverages
      % properties from both the `c_mat` and `c_mesh` classes and calls
      % `cal_cell_zgrad3d` to compute the gradient for each individual cell.
      %
      % Inputs:
      %   obj                   - (c_mat object) The instance of the c_mat class
      %                           containing material data and properties such as
      %                           material IDs, volume fractions, and mixed cell data.
      %   mesh                  - (c_mesh object) The instance of the c_mesh class
      %                           containing mesh-related properties such as cell
      %                           counts, boundary layers, and material information.
      %
      % Outputs:
      %   This function does not return any outputs but modifies the `vfgrad_mixcell`
      %   property of the `c_mat` class.
      %
      % Uses:
      %   obj.nmixcell             - (int32) The number of mixed cells in the problem.
      %   obj.nmat_in_mixcell      - (int array) The number of materials in each mixed zone.
      %   obj.matids_in_mixcell    - (cell array of int arrays) The material IDs in each
      %                              mixed zone.
      %   obj.vf_in_mixcell        - (cell array of double arrays) The volume fractions
      %                              of each material in each mixed zone.
      %   obj.ijk_in_mixcell       - (array of int arrays) The indices of each mixed
      %                              zone in the grid.
      %   mesh.ncell_prob          - (int array) The number of cells in the x, y,
      %                              and z directions, size: [3, 1].
      %   mesh.nbdry_prob          - (int) The number of boundary cells.
      %   mesh.xl_prob             - (double array) The lower bounds of the problem
      %                              domain in the x, y, and z directions, size: [3, 1].
      %   mesh.xr_prob             - (double array) The upper bounds of the problem
      %                              domain in the x, y, and z directions, size: [3, 1].
      %   mesh.nmat_for_3dcell     - (cell array of int arrays) The number of materials
      %                              in each 3D cell.
      %   mesh.matid_for_3dcell    - (cell array of int arrays) The material IDs in each
      %                              3D cell.
      %   mesh.mixcell_for_3dcell  - (cell array of int arrays) The mixed cell indices
      %                              for each 3D cell.
      %
      % Modifies:
      %   obj.vfgrad_mixcell       - (cell array of double arrays) The gradient of the
      %                              volume fractions for each material in each mixed zone.
      %
      % Original C declaration:
      % void cal_mixcell_zgrad3d(int *ncell, int nbdry, double *xl_prob, double *dx,
      %                          int ***nmat_for_3dcell, int ***matid_for_3dcell, int ***mixcell_for_3dcell,
      %                          int nmixzone, int *nmat_for_mzone,
      %                          int **matid_for_mzone, double **vf_for_mzone,
      %                          int **ijk_mzone, double ***vfgrad_for_mzone)

      % Initialize the last cell index accounting for boundary layers
      dim = 3;
      ncell_last = mesh.ncell_prob + 2 * mesh.nbdry_prob - 1;

      % Loop over all mixed zones
      for mix = 1:obj.nmixcell
        nm = obj.nmat_in_mixcell(mix);
        matids = obj.matids_in_mixcell{mix};
        vfs = obj.vf_in_mixcell{mix};
        ijk = obj.ijk_in_mixcell(mix,1:dim);

        for m = 1:nm
          i0 = ijk(1);
          j0 = ijk(2);
          k0 = ijk(3);

          % Skip boundary cells
          if i0 == 1 || i0 == ncell_last(1) ...
              || j0 == 1 || j0 == ncell_last(2) ...
              || k0 == 1 || k0 == ncell_last(3)
            continue;
          end

          % Initialize var and calculate the gradient
          matid = matids(m);
          var = zeros(3, 3, 3);
          var(2, 2, 2) = vfs(m);

          for k = -1:1
            kcell = k0 + k;
            k3 = k + 2;  % index in var[3][3][3]
            for j = -1:1
              jcell = j0 + j;
              j3 = j + 2;  % index in var[3][3][3]
              for i = -1:1
                if i == 0 && j == 0 && k == 0
                  continue;
                end
                icell = i0 + i;
                i3 = i + 2;  % index in var[3][3][3]

                var(k3, j3, i3) = 0.0;
                nm_cell = mesh.nmat_for_3dcell{kcell, jcell, icell};
                matid_cell = mesh.matid_for_3dcell{kcell, jcell, icell};
                if nm_cell == 1
                  if matid_cell == matid
                    var(k3, j3, i3) = 1.0;
                  end
                else
                  mix_nb = mesh.mixcell_for_3dcell{kcell, jcell, icell};
                  nm_nb = obj.nmat_in_mixcell(mix_nb);
                  assert(nm_nb == nm_cell);
                  matids_nb = obj.matids_in_mixcell{mix_nb};
                  vfs_nb = obj.vf_in_mixcell{mix_nb};
                  for m_nb = 1:nm_nb
                    if matid == matids_nb(m_nb)
                      var(k3, j3, i3) = vfs_nb(m_nb);
                      break;
                    end
                  end
                end
              end
            end
          end

          % Calculate the gradient for the current cell
          gradv = obj.cal_cell_zgrad3d(mesh.xr_prob - mesh.xl_prob, var);
          for i = 1:dim
            obj.vfgrad_mixcell{mix}(m, i) = gradv(i);
          end
        end
      end
    end

    function get_mpoly(obj, mesh)
      % get_mpoly - Generate material polygons or polyhedrons for mixed cells.
      %
      % This function generates material polygons (in 2D) or polyhedrons (in 3D)
      % for mixed cells within a simulation grid. The function calculates gradients,
      % determines the material interfaces, and reconstructs the geometry for each
      % mixed cell using properties from the `c_mat` and `c_mesh` classes.
      %
      % Inputs:
      %   obj                   - (c_mat object) The instance of the c_mat class
      %                           containing material data and properties such as
      %                           material IDs, volume fractions, and mixed cell data.
      %   mesh                  - (c_mesh object) The instance of the c_mesh class
      %                           containing mesh-related properties such as cell
      %                           counts, boundary layers, and material information.
      %
      % Outputs:
      %   This function does not return any outputs but modifies various properties
      %   of the `c_mat` class related to material polygons or polyhedrons.
      %
      % Uses:
      %   obj.nmixcell             - (int32) The number of mixed cells in the problem.
      %   obj.nmat_in_mixcell      - (int array) The number of materials in each mixed zone.
      %   obj.matids_in_mixcell    - (cell array of int arrays) The material IDs in each
      %                              mixed zone.
      %   obj.vf_in_mixcell        - (cell array of double arrays) The volume fractions
      %                              of each material in each mixed zone.
      %   obj.ijk_in_mixcell       - (array of int arrays) The indices of each mixed
      %                              zone in the grid.
      %   mesh.dim_prob            - (int) The dimensionality of the problem (2 or 3).
      %   mesh.ncell_prob          - (int array) The number of cells in the x, y,
      %                              and z directions, size: [3, 1] for 3D or [2, 1] for 2D.
      %   mesh.nbdry_prob          - (int) The number of boundary cells.
      %   mesh.xl_prob             - (double array) The lower bounds of the problem
      %                              domain in the x, y, and z directions, size: [3, 1] for 3D
      %                              or [2, 1] for 2D.
      %   mesh.xr_prob             - (double array) The upper bounds of the problem
      %                              domain in the x, y, and z directions, size: [3, 1] for 3D
      %                              or [2, 1] for 2D.
      %   mesh.nmat_for_2dcell     - (cell array of int arrays) The number of materials
      %                              in each 2D cell (used only if `dim_prob` is 2).
      %   mesh.matid_for_2dcell    - (cell array of int arrays) The material IDs in each
      %                              2D cell (used only if `dim_prob` is 2).
      %   mesh.mixcell_for_2dcell  - (cell array of int arrays) The mixed cell indices
      %                              for each 2D cell (used only if `dim_prob` is 2).
      %   mesh.nmat_for_3dcell     - (cell array of int arrays) The number of materials
      %                              in each 3D cell (used only if `dim_prob` is 3).
      %   mesh.matid_for_3dcell    - (cell array of int arrays) The material IDs in each
      %                              3D cell (used only if `dim_prob` is 3).
      %   mesh.mixcell_for_3dcell  - (cell array of int arrays) The mixed cell indices
      %                              for each 3D cell (used only if `dim_prob` is 3).
      %
      % Modifies:
      %   obj.mix_mpoly_to_mix     - (int array) Mapping of mixed cell indices to
      %                              multipolygon indices.
      %   obj.nnode_in_mixcell     - (int array) The number of nodes in each mixed cell
      %                              polygon or polyhedron.
      %   obj.coords_in_mixcell    - (cell array of double arrays) Coordinates of nodes
      %                              in each mixed cell.
      %   obj.nnode_for_minterface_in_mixcell - (cell array of int arrays) The number of nodes
      %                                         for each material interface in each mixed cell.
      %   obj.nodes_for_minterface_in_mixcell - (cell array of int arrays) The nodes
      %                                         for each material interface in each mixed cell.
      %   obj.nnode_for_mpoly_in_mixcell      - (cell array of int arrays, 2D only) The number of
      %                                         nodes for each material polygon in each mixed cell.
      %   obj.nodes_for_mpoly_in_mixcell      - (cell array of int arrays, 2D only) The nodes
      %                                         for each material polygon in each mixed cell.
      %   obj.nface_for_mpoly_in_mixcell      - (cell array of int arrays, 3D only) The number of
      %                                         faces for each material polyhedron in each mixed cell.
      %   obj.nnode_for_face_ea_mpoly_in_mixcell - (cell array of int arrays, 3D only) The number of
      %                                            nodes for each face of a material polyhedron.
      %   obj.nodelist_for_face_ea_mpoly_in_mixcell - (cell array of int arrays, 3D only) The nodes
      %                                               for each face of a material polyhedron.
      %
      % Original C declaration:
      % void get_mpoly(int dim, int *ncell, int nbdry, double *xl_prob, double *dx,
      %                int  **nmat_for_2dcell, int  **matid_for_2dcell, int  **mixcell_for_2dcell,
      %                int ***nmat_for_3dcell, int ***matid_for_3dcell, int ***mixcell_for_3dcell)

      global vof2d;

      % Initialize variables
      dim = mesh.dim_prob;
      geop = 1; % Cartesian coordinate
      nmixcell = obj.nmixcell;

      if nmixcell == 0
        return;
      end

      % Initialize cell arrays
      ncell_ext = mesh.ncell_prob + 2 * mesh.nbdry_prob;
      ncell_last = ncell_ext(1:dim);

      % Set up storage
      obj.set_mpoly_storage(dim);

      % Calculate volume
      vol = prod(mesh.xr_prob - mesh.xl_prob);

      % Allocate memory for vfgrad
      vfgrad_for_mixcell = cell(nmixcell, 1);
      nmsum = sum(obj.nmat_in_mixcell);
      nmmax = max(obj.nmat_in_mixcell);

      % Assign vfgrad cells
      for mix = 1:nmixcell
        nm = obj.nmat_in_mixcell(mix);
        vfgrad_for_mixcell{mix} = zeros(nm, dim);
      end

      % Calculate gradients
      if dim == 2
        obj.cal_mixcell_zgrad2d(mesh);
      elseif dim == 3
        obj.cal_mixcell_zgrad3d(mesh);
      end

      % Allocate memory for normals and priorities
      normal_ea_mat = zeros(nmmax, dim);
      priority_ea_mat = zeros(nmmax, 1);

      obj.mix_mpoly_to_mix = zeros(nmixcell, 1);
      obj.mix_to_mpoly_mix = zeros(nmixcell, 1);

      if dim == 2
        sizes = [nmmax, 3, 3];
        nmat_2dsmesh = zeros(3, 3);
        matid_2dsmesh = cell(3, 3);
        vf_2dsmesh = cell(3, 3);

      elseif dim == 3
        sizes = [3, 3, 3];
        nmat_3dsmesh = zeros(3, 3, 3);
        matid_3dsmesh = cell(3, 3, 3);
        vf_3dsmesh = cell(3, 3, 3);
      end

      nmixcell_mpoly = 1;

      for mix = 1:nmixcell
        nm = obj.nmat_in_mixcell(mix);
        ijk = obj.ijk_in_mixcell(mix,1:dim);

        % Skip boundary cells
        if any(ijk == 1 | ijk == ncell_last)
          continue;
        end

        xl = mesh.xl_prob + double(ijk - mesh.nbdry_prob - 1) ...
          .* (mesh.xr_prob - mesh.xl_prob) ...
          ./ (mesh.ncell_prob);

        % Determine normal vectors and priorities
        vfs = obj.vf_in_mixcell{mix};
        for m = 1:nm
          gradvf = vfgrad_for_mixcell{mix}(m,:);
          norm = -gradvf;
          vsq = sum(norm.^2);
          norm = norm / sqrt(vsq);
          priority_ea_mat(m) = sqrt(vsq) * sqrt(vfs(m));
        end

        [~, m0] = max(priority_ea_mat);
        norm0 = normal_ea_mat(m0, :);

        % Invert normals if necessary
        if mod(m0, 2) == 1
          normal_ea_mat = -norm0;
        end

        nnode_tot = 0;
        coords_tot = [];

        if dim == 2
          for j0 = 1:3
            j = ijk(2) + j0 - 2;
            for i0 = 1:3
              i = ijk(1) + i0 - 2;
              mynm = mesh.nmat_for_2dcell(j, i);
              nmat_2dsmesh(j0, i0) = mynm;
              if mynm == 1
                matid_2dsmesh{j0, i0} = mesh.matid_for_2dcell(j, i);
                vf_2dsmesh{j0, i0} = 1.0;
              else
                mymix = mesh.mixcell_for_2dcell(j, i);
                matid_2dsmesh{j0, i0} = obj.matids_in_mixcell{mymix};
                vf_2dsmesh{j0, i0} = obj.vf_in_mixcell{mymix};
              end
            end
          end

          [coords_tot, nnode_tot, nodes_for_minterface, ...
            nodes_for_mpoly, nnode_for_minterface, nnode_for_mpoly] ...
            = vof2d.reconstruct2d_nmat_pagosa(geop, xl, mesh.dx, ...
            nmat_2dsmesh, matid_2dsmesh, vf_2dsmesh);

        elseif dim == 3
          for k0 = 1:3
            k = ijk(3) + k0 - 2;
            for j0 = 1:3
              j = ijk(2) + j0 - 2;
              for i0 = 1:3
                i = ijk(1) + i0 - 2;
                mynm = mesh.nmat_for_3dcell{k, j, i};
                nmat_3dsmesh(k0, j0, i0) = mynm;
                if mynm == 1
                  matid_3dsmesh{k0, j0, i0} = mesh.matid_for_3dcell{k, j, i};
                  vf_3dsmesh{k0, j0, i0} = 1.0;
                else
                  mymix = mesh.mixcell_for_3dcell{k, j, i};
                  matid_3dsmesh{k0, j0, i0} = obj.matids_in_mixcell{mymix};
                  vf_3dsmesh{k0, j0, i0} = obj.vf_in_mixcell{mymix};
                end
              end
            end
          end

          [coords_tot, nnode_tot, nnode_for_minterface, ...
            nodes_for_minterface, ...
            nface_for_mpoly, nnode_for_face_ea_mpoly, nodelist_for_face_ea_mpoly] = ...
            obj.reconstruct3d_nmat_pagosa(xl, mesh.xr_prob - mesh.xl_prob, ...
            nmat_3dsmesh, matid_3dsmesh, vf_3dsmesh);
        end

        % Store results
        obj.mix_mpoly_to_mix(nmixcell_mpoly) = mix;
        obj.nnode_in_mixcell(nmixcell_mpoly) = nnode_tot;
        obj.coords_in_mixcell{nmixcell_mpoly} = coords_tot;
        obj.nnode_for_minterface_in_mixcell{nmixcell_mpoly} = nnode_for_minterface;
        obj.nodes_for_minterface_in_mixcell{nmixcell_mpoly} = nodes_for_minterface;

        if dim == 2
          obj.nnode_for_mpoly_in_mixcell{nmixcell_mpoly} = nnode_for_mpoly;
          obj.nodes_for_mpoly_in_mixcell{nmixcell_mpoly} = nodes_for_mpoly;
        elseif dim == 3
          obj.nface_for_mpoly_in_mixcell{nmixcell_mpoly} = nface_for_mpoly;
          obj.nnode_for_face_ea_mpoly_in_mixcell{nmixcell_mpoly} = nnode_for_face_ea_mpoly;
          obj.nodelist_for_face_ea_mpoly_in_mixcell{nmixcell_mpoly} = nodelist_for_face_ea_mpoly;
        end

        obj.mix_to_mpoly_mix(mix) = nmixcell_mpoly;
        nmixcell_mpoly = nmixcell_mpoly + 1;
      end
      obj.nmixcell_mpoly = nmixcell_mpoly - 1;

      % Free allocated memory
      if dim == 2
        clear nmat_2dsmesh matid_2dsmesh vf_2dsmesh;
      elseif dim == 3
        clear nmat_3dsmesh matid_3dsmesh vf_3dsmesh;
      end
    end

    function [nmix_cell] = set_mat(obj, mesh, btype_lower, btype_upper, ...
        is_solid, gamma_ea_mat, nreg, reg2matids, reg_shape, rho_ea_reg, ...
        pres_ea_reg, ei_ea_reg, v_ea_reg)
      % Function: set_mat
      %
      % This function sets up the material properties for a mesh, depending on
      % the dimensionality (2D or 3D). It processes boundary conditions and
      % determines the properties for each cell in the mesh.
      %
      % Inputs:
      % - obj (c_mat): The current c_mat object instance
      % - mesh (c_mesh): An instance of the c_mesh class
      % - btype_lower (bdry_type array): Boundary type at the lower end of each
      %   dimension
      % - btype_upper (bdry_type array): Boundary type at the upper end of each
      %   dimension
      % - is_solid (int array): Indicator array for solid materials
      % - gamma_ea_mat (double array): Gamma values for each material
      % - nreg (int): Number of regions
      % - reg2matids (int array): Material IDs for each region
      % - reg_shape (region_shape array): Geometric shapes defining regions
      % - rho_ea_reg (double array): Density for each region
      % - pres_ea_reg (double array): Pressure for each region
      % - ei_ea_reg (double array): Internal energy density for each region
      % - v_ea_reg (cell): Cell array containing velocity arrays for each region
      %
      % Outputs:
      % - nmix_cell (int): Number of mixed cells
      %
      % Uses:
      % - mesh.dim_prob (int): Dimensionality of the problem (2D or 3D)
      % - mesh.dx (double array): Cell sizes in each dimension
      % - mesh.ncell_prob (int array): Number of cells in each dimension
      % - mesh.xl_prob (double array): Lower bounds of the simulation domain
      % - mesh.nbdry_prob (int): Number of boundary layers
      % - mesh.nmat_mesh (int): Number of materials in the mesh
      % - mesh.matids_mesh (int array): Material IDs in the mesh
      %
      % Modifies:
      % - obj.nmat_prob (int): Number of materials in the problem
      % - obj.matids_prob (int array): Material IDs for each material in the
      %   problem
      % - mesh.nmat_for_2dcell (int array): Number of materials in each 2D cell
      % - mesh.matid_for_2dcell (int array): Material IDs for each 2D cell
      % - mesh.mixcell_for_2dcell (int array): Mixed cell indices for 2D cells
      % - mesh.rho_for_2dcell (double array): Density for each 2D cell
      % - mesh.ei_for_2dcell (double array): Internal energy density for each
      %   2D cell
      % - mesh.pres_for_2dcell (double array): Pressure for each 2D cell
      % - mesh.vel_for_2dcell (double array): Velocity for each 2D cell
      % - mesh.nmat_for_3dcell (int array): Number of materials in each 3D cell
      % - mesh.matid_for_3dcell (int array): Material IDs for each 3D cell
      % - mesh.mixcell_for_3dcell (int array): Mixed cell indices for 3D cells
      % - mesh.rho_for_3dcell (double array): Density for each 3D cell
      % - mesh.ei_for_3dcell (double array): Internal energy density for each
      %   3D cell
      % - mesh.pres_for_3dcell (double array): Pressure for each 3D cell
      % - mesh.vel_for_3dcell (double array): Velocity for each 3D cell
      %
      % Original C declaration:
      % void set_mat(int dim, int *ncell, double *xl_prob, double *dx, int nbdry,
      %              Bdry_Type *btype_lower, Bdry_Type *btype_upper,
      %              int nmat, int *matids_ea_mat, int *is_solid, double *gamma_ea_mat,
      %              int nreg, int *reg2matids, Region_Shape *reg_shape,
      %              double *rho_ea_reg, double *pres_ea_reg, double *ei_ea_reg, double **v_ea_reg,
      %              int *nmix_cell,
      %              int **nmat_for_2dcell, int **matid_for_2dcell, int **mixcell_for_2dcell,
      %              double **rho_for_2dcell, double **ei_for_2dcell, double **pres_for_2dcell, double ***vel_for_2dcell,
      %              int ***nmat_for_3dcell, int ***matid_for_3dcell, int ***mixcell_for_3dcell,
      %              double ***rho_for_3dcell, double ***ei_for_3dcell, double ***pres_for_3dcell, double ****vel_for_3dcell)

      global bdry;
      dim = mesh.dim_prob;
      dx = mesh.dx;

      obj.nmat_prob = mesh.nmat_mesh;
      obj.matids_prob = mesh.matids_mesh;

      if dim == 2
        obj.set_2dmesh_mat(mesh, nreg, is_solid, gamma_ea_mat, ...
          reg2matids, reg_shape, rho_ea_reg, pres_ea_reg, ...
          ei_ea_reg, v_ea_reg);

      elseif dim == 3
        obj.set_3dmesh_mat(mesh, nreg, is_solid, gamma_ea_mat, ...
          reg2matids, reg_shape, rho_ea_reg, pres_ea_reg, ...
          ei_ea_reg, v_ea_reg);
      end

      % Apply boundary conditions to mixed cells
      bdry.mixcell_bdry_condition(mesh, obj, btype_lower, btype_upper);


      % Apply boundary conditions to cell variables for velocity
      if dim == 2
        bdry.bdry_cell_2d(mesh, btype_lower, btype_upper);
        bdry.bdry_cell_vel_2d(mesh, btype_lower, btype_upper);
      elseif dim == 3
        bdry.bdry_cell_3d(mesh, btype_lower, btype_upper);
        bdry.bdry_cell_vel_3d(mesh, btype_lower, btype_upper);
      end

      nmix_cell = obj.nmixcell;  % Return the updated nmix_cell

    end

    function mat_sound_speed(obj, mesh, nmat, matids, is_solid, gamma_ea_mat)
      % mat_sound_speed Calculates the sound speed for mixed cells in the mesh.
      %
      % This method computes the sound speed for each 2D or 3D mixed cell in the mesh
      % based on the materials present in the cells. It handles both solid and non-solid
      % materials, applying different methods to compute the sound speed.
      %
      % Inputs:
      %   obj          - (c_mat) The current material object instance.
      %   mesh         - (c_mesh) The mesh object containing grid and velocity data.
      %   nmat         - (int) Number of materials in the mesh.
      %   matids       - (array of int) Material IDs for each material.
      %   is_solid     - (array of int) Indicator array for solid materials.
      %   gamma_ea_mat - (array of double) Gamma values for each material.
      %
      % Outputs:
      %   None. The function updates the sound speed arrays for 2D and 3D cells.
      %
      % Uses:
      %   mesh.ncell_prob         - (array of int) Number of cells in each spatial dimension.
      %   mesh.nbdry_prob         - (int) Number of boundary cells.
      %   obj.nmat_in_mixcell     - (array of int) Number of materials in each mixed cell.
      %   obj.matids_in_mixcell   - (cell of int arrays) Material IDs in each mixed cell.
      %   obj.rho_in_mixcell      - (cell of double arrays) Density for each material in mixed cells.
      %   obj.pres_in_mixcell     - (cell of double arrays) Pressure for each material in mixed cells.
      %   obj.vf_in_mixcell       - (cell of double arrays) Volume fraction for each material in mixed cells.
      %   obj.ijk_in_mixcell      - (cell of int arrays) Indices of mixed cells in the mesh.
      %
      % Modifies:
      %   mesh.cs_for_2dcell      - (2D array of double) Updates sound speed for 2D cells.
      %   mesh.cs_for_3dcell      - (3D array of double) Updates sound speed for 3D cells.
      %
      % Original C declaration:
      % void mat_sound_speed(int dim, int *ncell, int nbdry,
      %                      int nmat, int *matids, int *is_solid, double *gamma_ea_mat,
      %                      double **cs_for_2dcell,  double ***cs_for_3dcell);

      global eos tiny;  % Global object for EOS (Equation of State)

      dim = mesh.dim_prob;       % Dimensionality of the mesh
      ncell = mesh.ncell_prob;   % Number of cells in each direction
      nbdry = mesh.nbdry_prob;   % Number of boundary cells

      % Initialize extended cell size
      ncell_ext = ncell + 2 * nbdry;

      matid = max(obj.matids_prob(1:obj.nmat_prob));
      gamma_matid = zeros(matid, 1);

      for i = 1:nmat
        gamma_matid(matids(i)) = gamma_ea_mat(i);
      end

      % Initialize the sound speed array for mixed cells
      cs = zeros(obj.nmixcell, 1);

      % Calculate sound speeds for each mixed cell
      for mx = 1:obj.nmixcell
        nm = obj.nmat_in_mixcell(mx);
        for m = 1:nm
          matid = obj.matids_in_mixcell{mx}(m);
          if ~is_solid(matid)
            csm = sqrt(gamma_matid(matid) * obj.pres_in_mixcell{mx}(m) / ...
              (obj.rho_in_mixcell{mx}(m) + tiny));
          else
            csm = eos.sound_speed_solid(obj.rho_in_mixcell{mx}(m));
          end
          cs(mx) = cs(mx) + obj.vf_in_mixcell{mx}(m) * csm;
        end
      end

      % Assign sound speeds to cells (2D and 3D cases)
      if dim == 2
        for mx = 1:obj.nmixcell
          ijk = obj.ijk_in_mixcell(mx,1:dim);
          i = ijk(1);
          j = ijk(2);
          mesh.cs_for_2dcell(j, i) = cs(mx);
        end
      elseif dim == 3
        for mx = 1:obj.nmixcell
          ijk = obj.ijk_in_mixcell(mx,1:dim);
          i = ijk(1);
          j = ijk(2);
          k = ijk(3);
          mesh.cs_for_3dcell(k, j, i) = cs(mx);
        end
      end
    end

    function [nmat_advected, matid_advected, vol_advected, mass_advected, ener_advected] = advect2d(...
        obj, xl_cell, dx, mixcell, xl_slab, xr_slab, inward_norm, plane_of_slab)
      % advect2d Performs material advection across a slab in a 2D mixed cell.
      %
      % This method calculates the material volume, mass, and energy advected
      % across a slab interface within a 2D mixed cell, based on the material's
      % polygon geometry and the slab's orientation.
      %
      % Inputs:
      %   obj             - (c_mat) The current material object instance.
      %   xl_cell         - (array of double) Coordinates of the lower corner of the cell [size: 1x2].
      %   dx              - (array of double) Cell dimensions [size: 1x2].
      %   mixcell         - (int) Index of the mixed cell.
      %   xl_slab         - (array of double) Coordinates of the lower corner of the slab [size: 1x2].
      %   xr_slab         - (array of double) Coordinates of the upper corner of the slab [size: 1x2].
      %   inward_norm     - (array of double) Inward normal vector of the slab [size: 1x2].
      %   plane_of_slab   - (int) Identifier for the slab plane (0 for xl, 1 for xr, 2 for yl, 3 for yr).
      %
      % Outputs:
      %   nmat_advected   - (int) Number of materials advected across the slab.
      %   matid_advected  - (array of int) Material IDs for the advected materials.
      %   vol_advected    - (array of double) Volumes of the advected materials.
      %   mass_advected   - (array of double) Masses of the advected materials.
      %   ener_advected   - (array of double) Energies of the advected materials.
      %
      % Uses:
      %   obj.mix_to_mpoly_mix   - (array of int) Mapping of mixed cells to material polygons.
      %   obj.coords_in_mixcell  - (cell of arrays) Coordinates of material polygons in each mixed cell.
      %   obj.nnode_in_mixcell   - (array of int) Number of nodes in each mixed cell.
      %   obj.nmat_in_mixcell    - (array of int) Number of materials in each mixed cell.
      %   obj.matids_in_mixcell  - (cell of arrays) Material IDs in each mixed cell.
      %   obj.rho_in_mixcell     - (cell of arrays) Densities for each material in mixed cells.
      %   obj.ei_in_mixcell      - (cell of arrays) Internal energies for each material in mixed cells.
      %
      % Modifies:
      %   None.
      %
      % Original C declaration:
      % void advect2d(double *xl_cell, double *dx_cell, int mixcell,
      %               double *xl_slab, double *xr_slab, double *inward_norm, int plane_of_slab,
      %               int *nmat_advected, int *matid_advected, double *vol_advected,
      %               double *mass_advected, double *ener_advected)
      global vof2d util;
      dim = 2; dx = dx(1:dim);
      accuracy = 1.0e-06; % accuracy used to check below or above the interface
      factor_vol = prod(dx);  % Volume scaling factor

      % Scale slab coordinates
      xl_slab_scaled = (xl_slab - xl_cell) ./ dx;
      xr_slab_scaled = (xr_slab - xl_cell) ./ dx;

      % Set interface coordinates based on plane of slab
      if plane_of_slab == 0  % xl-plane of the slab
        c0 = [xl_slab_scaled(1), xr_slab_scaled(2)];
        c1 = [xl_slab_scaled(1), xl_slab_scaled(2)];
      elseif plane_of_slab == 1  % xr-plane of the slab
        c0 = [xr_slab_scaled(1), xl_slab_scaled(2)];
        c1 = [xr_slab_scaled(1), xr_slab_scaled(2)];
      elseif plane_of_slab == 2  % yl-plane of the slab
        c0 = [xl_slab_scaled(1), xl_slab_scaled(2)];
        c1 = [xr_slab_scaled(1), xl_slab_scaled(2)];
      elseif plane_of_slab == 3  % yr-plane of the slab
        c0 = [xr_slab_scaled(1), xr_slab_scaled(2)];
        c1 = [xl_slab_scaled(1), xr_slab_scaled(2)];
      end

      % Calculate distance from the slab
      distance = 0.5*(inward_norm(1)*(c0(1) + c1(1)) + inward_norm(2)*(c0(2) + c1(2)));

      % Retrieve the mixed cell polygon
      mx_mpoly = obj.mix_to_mpoly_mix(mixcell);
      assert(mx_mpoly >= 1 && mx_mpoly <= obj.nmixcell_mpoly);

      nnode = obj.nnode_in_mixcell(mx_mpoly);
      coords_intersect = obj.coords_in_mixcell{mx_mpoly};  % Coordinates of material polygon
      coords_intersect = reshape(coords_intersect, [2, nnode])';

      % Scale the coordinates of the material polygons
      coords_scaled = coords_intersect;
      for n = 1:nnode
        coords_scaled(n, :) = (coords_intersect(n, :) - xl_cell) ./ dx;
      end

      % Initialize variables for advection
      nmat_advected = 0;
      vol_advected = [];
      mass_advected = [];
      ener_advected = [];
      matid_advected = [];

      nm = obj.nmat_in_mixcell(mixcell);
      assert(nm > 1);

      for m = 1:nm
        nn = obj.nnode_for_mpoly_in_mixcell{mx_mpoly}(m);
        nodes = obj.nodes_for_mpoly_in_mixcell{mx_mpoly}{m};

        % Copy the coordinates of the material polygon to coords_mpoly
        coords_mpoly = coords_scaled(nodes, :);

        % Compute distances of nodes from the interface
        ds_ea_node = inward_norm(1) * coords_mpoly(:, 1) + inward_norm(2) * coords_mpoly(:, 2);

        % Determine node locations relative to the interface
        node_loc = zeros(nn, 1);
        node_loc(abs(ds_ea_node - distance) <= accuracy) = 0;  % on the plane
        node_loc(ds_ea_node > distance) = 1;  % above the plane
        node_loc(ds_ea_node < distance) = -1;  % below the plane

        % Find the interface and split the material polygon
        coords_mpoly_flat = coords_mpoly';
        coords_mpoly_flat = coords_mpoly_flat(:);
        [~, coords_new, ~, ~, ~, nodelist_lower, nnode_upper, nodelist_upper] = ...
          vof2d.find_interface2d(nn, coords_mpoly_flat, (1:4), inward_norm, distance, node_loc);
        coords_new = reshape(coords_new, 2, [])';


        % If upper portion has more than 2 nodes, calculate volume, mass, and energy
        if nnode_upper > 2
          coords_intersect = zeros(2, nnode_upper);
          for idx = 1:nnode_upper
            n = nodelist_upper(idx);
            if (n <= nn)
              coords_intersect(:, idx) = coords_mpoly(n, :);
            else
              coords_intersect(:, idx) = coords_new(n - nn, :);
            end
          end
          vol = util.cal_poly_area(nnode_upper, coords_intersect(:), nnode_upper, []) * factor_vol;

          % Store advected material data
          vol_advected(end + 1) = vol;
          matid_advected(end + 1) = obj.matids_in_mixcell{mixcell}(m);
          mass_advected(end + 1) = obj.rho_in_mixcell{mixcell}(m) * vol;
          ener_advected(end + 1) = obj.ei_in_mixcell{mixcell}(m) * vol;

          nmat_advected = nmat_advected + 1;
        end

        if ~isempty(nodelist_lower), nodelist_lower = []; end
        if ~isempty(nodelist_upper), nodelist_upper = []; end
      end
    end

    function [nmat_int, matids_int, matidx_int, vol_int] = intersect2d(obj, xl_cell, dx, mixcell)
      % intersect2d Finds the intersection of a cell with material polygons in a mixed cell.
      %
      % This function calculates the intersection of a rectangular cell defined by
      % `xl_cell` and `dx` with the material polygons in a mixed cell, determining
      % the number of materials, their IDs, and the corresponding volumes.
      %
      % Inputs:
      %   obj       - (c_bdry) The current boundary object instance.
      %   xl_cell   - (array of double) Coordinates of the lower corner of the cell [size: 1x2].
      %   dx        - (array of double) Cell dimensions [size: 1x2].
      %   mixcell   - (int) Index of the mixed cell.
      %
      % Outputs:
      %   nmat_int  - (int) Number of intersecting materials.
      %   matids_int - (array of int) Material IDs for intersecting materials.
      %   matidx_int - (array of int) Indices of the intersecting materials in the mixed cell.
      %   vol_int   - (array of double) Volumes of the intersecting materials.
      %
      % Uses:
      %   obj.mix_to_mpoly_mix      - (array of int) Mapping of mixed cells to material polygons.
      %   obj.nmat_in_mixcell       - (array of int) Number of materials in each mixed cell.
      %   obj.matids_in_mixcell     - (cell array of int arrays) Material IDs for each mixed cell.
      %   obj.nnode_for_mpoly_in_mixcell - (cell array of int arrays) Number of nodes for each material polygon.
      %   obj.nodes_for_mpoly_in_mixcell - (cell array of int arrays) Node indices for each material polygon.
      %   obj.coords_in_mixcell     - (cell array of double arrays) Coordinates of material polygons in each mixed cell.
      %
      % Modifies:
      %   None.
      %
      % Original C declaration:
      % void intersect2d(double *xl_cell, double *dx, int mixcell,
      %                  int *nmat_int, int *matids_int, int *matidx_int, double *vol_int)

      dim = 2;
      vol_cell = prod(dx);  % Volume of the cell
      nnode_rec = 4;

      % Create rectangle coordinates (cell)
      coords_rec = [
        xl_cell(1), xl_cell(2);
        xl_cell(1) + dx(1), xl_cell(2);
        xl_cell(1) + dx(1), xl_cell(2) + dx(2);
        xl_cell(1), xl_cell(2) + dx(2)
        ];

      % Get the mixed cell polygon
      mpoly_mix = obj.mix_to_mpoly_mix(mixcell);
      assert(mpoly_mix >= 0);

      nmat = obj.nmat_in_mixcell(mixcell);
      matids = obj.matids_in_mixcell{mixcell};
      nnode_for_mpoly = obj.nnode_for_mpoly_in_mixcell{mpoly_mix};
      nodes_for_mpoly = obj.nodes_for_mpoly_in_mixcell{mpoly_mix};
      coords_one_mixcell = obj.coords_in_mixcell{mpoly_mix};

      % Initialize output variables
      nmat_int = 0;
      matids_int = [];
      matidx_int = [];
      vol_int = [];

      % Iterate through each material polygon in the mixed cell
      for idx = 1:nmat
        nn_mpoly = nnode_for_mpoly(idx);
        nodes = nodes_for_mpoly{idx};

        % Copy material polygon coordinates
        coords_mpoly = zeros(nn_mpoly, dim);
        for k = 1:nn_mpoly
          n = nodes(k);
          coords_mpoly(k, :) = coords_one_mixcell(n, :);
        end

        % Find intersection with the rectangular cell
        [nnode_int, coords_int] = remap2d_scaled(0, nnode_rec, coords_rec, ...
          inward_norm_ea_face_default, nn_mpoly, coords_mpoly);

        % Calculate the volume of the intersecting polygon
        if nnode_int > 2
          vol = cal_poly_area(nnode_int, coords_int) * vol_cell;
          if vol / vol_cell > small
            nmat_int = nmat_int + 1;
            matids_int(nmat_int) = matids(idx);
            matidx_int(nmat_int) = idx;
            vol_int(nmat_int) = vol;
          end
        end
      end
    end



    %%%%%%%%%%%%%%% empty stubs
    function mat_pass_mix(obj, nmix, ijk_for_mixcell, nmat_for_mixcell, matids_for_mixcell, ...
        nnode_for_mixcell, coords_for_mixcell, ...
        nnode_for_mpoly_for_mixcell, nodes_for_mpoly_for_mixcell, ...
        nface_for_mpoly_for_mixcell, nnode_for_face_ea_mpoly_for_mixcell, ...
        nodelist_for_face_ea_mpoly_for_mixcell)
      % TODO: Implement mat_pass_mix
    end

    function mat_sort_mixcells(obj, dim, ncell, nbdry)
      % TODO: Implement mat_sort_mixcells
    end

    function [nummixcell, nummixcell_int, nummixcell_mpoly, nmat_ea_mixcell, ...
        ijk_ea_mixcell, matids_ea_mixcell, vf_ea_mixcell, ...
        rho_ea_mixcell, pres_ea_mixcell, ei_ea_mixcell] ...
        = mat_pass_mix_data(obj)
      % TODO: Implement mat_pass_mix_data
      % OK: not sure this function is necessary.
    end

    function mat_set_mix_data(obj, nmixcell_new, nmixcell_int_new, ...
        nmat_in_mixcell_new, ijk_in_mixcell_new, matids_in_mixcell_new, ...
        vf_in_mixcell_new, rho_in_mixcell_new, ...
        pres_in_mixcell_new, ei_in_mixcell_new)
      % TODO: Implement mat_set_mix_data
    end

    function mat_write_mat(obj, fileid, dim, pass)
      % TODO: Implement mat_write_mat
    end

    function clean_mpoly(obj)
      % Dummy implementation of clean_mpoly for completeness
      % In the actual class, this method should clean up existing storage.
    end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%% debug helpers
    function [nmat_cell, matids_cell, rho, pres, ei, vel, ...
        vf_ea_mat, rho_ea_mat, pres_ea_mat, ei_ea_mat] = ...
        debug_set_cell_mat(obj, dim, xl, dx, nreg, is_solid, ...
        matids_ea_reg, gamma_ea_mat, reg_shape, ...
        rho_ea_reg, pres_ea_reg, ei_ea_reg, v_ea_reg)
      % DEBUG_SET_CELL_MAT Debugs the inputs, outputs, and uses of set_cell_mat.
      %
      % This function prints all the input, output, and internal variables
      % before and after calling the set_cell_mat function. It also prints
      % the details of the region_shape class.

      % Display inputs before the function call
      disp('Before calling set_cell_mat:');

      % Display scalars and arrays
      fprintf('dim = %d\n', dim);
      fprintf('xl = [%f, %f]\n', xl(1), xl(2));
      fprintf('dx = [%f, %f]\n', dx(1), dx(2));

      % Material-related inputs
      fprintf('obj.nmat = %d\n', obj.nmat_prob);
      fprintf('matids_ea_mat = %s\n', mat2str(obj.matids_prob));
      fprintf('is_solid = %s\n', mat2str(is_solid));
      fprintf('gamma_ea_mat = %s\n', mat2str(gamma_ea_mat));
      fprintf('nreg = %d\n', nreg);
      fprintf('matids_ea_reg = %s\n', mat2str(matids_ea_reg));

      % Region properties (rho, pressure, energy)
      fprintf('rho_ea_reg = %s\n', mat2str(rho_ea_reg));
      fprintf('pres_ea_reg = %s\n', mat2str(pres_ea_reg));
      fprintf('ei_ea_reg = %s\n', mat2str(ei_ea_reg));

      % Velocities for each region
      disp('v_ea_reg:');
      for i = 1:nreg
        fprintf('  v_ea_reg{%d} = [%f, %f, %f]\n', i, v_ea_reg{i}(1), v_ea_reg{i}(2), v_ea_reg{i}(3));
      end

      % Print region_shape class properties
      disp('reg_shape:');
      for i = 1:nreg
        fprintf('  reg_shape{%d}.type = %s\n', i, char(reg_shape(i).type));  % Assuming type is enum Shape_Type
        fprintf('  reg_shape{%d}.parameters = %s\n', i, mat2str(reg_shape(i).parameters));
      end

      % Call the actual set_cell_mat function and capture outputs
      [nmat_cell, matids_cell, rho, pres, ei, vel, ...
        vf_ea_mat, rho_ea_mat, pres_ea_mat, ei_ea_mat] = ...
        set_cell_mat(obj, dim, xl, dx, nreg, is_solid, ...
        matids_ea_reg, gamma_ea_mat, reg_shape, ...
        rho_ea_reg, pres_ea_reg, ei_ea_reg, v_ea_reg);

      % Display outputs after the function call
      disp('After calling set_cell_mat:');

      % Output results
      fprintf('nmat_cell = %d\n', nmat_cell);
      fprintf('matids_cell = %s\n', mat2str(matids_cell));

      fprintf('vf_ea_mat = %s\n', mat2str(vf_ea_mat));
      fprintf('rho_ea_mat = %s\n', mat2str(rho_ea_mat));
      fprintf('pres_ea_mat = %s\n', mat2str(pres_ea_mat));
      fprintf('ei_ea_mat = %s\n', mat2str(ei_ea_mat));

      fprintf('rho = %f\n', rho);
      fprintf('pres = %f\n', pres);
      fprintf('ei = %f\n', ei);
      fprintf('vel = [%f, %f]\n', vel(1), vel(2));

    end

  end
end
