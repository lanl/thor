classdef SInit2 < handle
  %INIT_STRUCT Summary of this class goes here
  %   Detailed explanation goes here

  properties
    filename = '';        % name of the parameter file
    dims = 0;             % number of problem dimensions
    ncells = [0, 0, 0];   % array(1:dims): mesh sizes
    nbdry = 0;            % (small positive integer): boundary zone size
    btype_lower bdry_type % lower boundary
    btype_upper bdry_type % upper boundary
    xl_prob = [0, 0, 0];  % [1:dims]: left-bottom corner of the domain
    xr_prob = [0, 0, 0];  % [1:dims]: top-right corner of the domain
    nummat = 1;           % number of materials
    matids_ea_mat = [0];  % [1:nummat]: material IDs for each material
    is_solid = [0];       % [1:nummat]: is each material solid or not
    gamma_ea_mat = [0];   % [1:nummat]: polytropic gamma for each material
    tmax = 0.0;           % final time
    ncycle = 0;           % cycle number (iteration)
    dt_viz = 0.0;         % how often to dump visualization data
    ncycle_viz = 1;       % how often to dump viz. data in cycles
    courant = 0.5;        % Courant factor limit
    dt = 0.1;             % timestep size
    problemname = '';     % problem name
  end
  methods
    function [obj, status] = SInit2(filename)
      % Constructor method

      % Constants
      name_accelpiston = 'accelpiston';
      name_shockformation = 'shockformation';
      name_blowoff = 'blowoff';
      name_advection = 'advection';
      name_2shocks = '2shocks';
      name_sod = 'sod';
      name_mach60 = 'mach60';
      nmat_mx = 4; % max number of materials
      nreg_mx = 6; % max number of regions

      % Initialize variables
      % RC: Towards removing clutter commented out are variables
      % that do not need preallocation of storage
      % inputtype = '';
      % equalsigne = '';
      % probname = '';
      % dim = 0;
      % nmat = 0;
      nreg = 0;
      nb = 2;
      dir = 0;
      reg2matids = zeros(1, nreg_mx);
      ncell = zeros(1, 3);
      xl = zeros(1, 3);
      xr = zeros(1, 3);
      rho_ea_reg = zeros(1, nreg_mx);
      pres_ea_reg = zeros(1, nreg_mx);
      ei_ea_reg = zeros(1, nreg_mx);
      v_ea_reg = zeros(nreg_mx, 3);
      reg_shape = repmat(region_shape(shape_type.shape_universe, []), 1, nreg_mx);  % Array of region_shape objects
      vel_ea_reg = cell(1, nreg_mx);
      param_left = zeros(1, 6);
      param_block = zeros(1, 6);
      param_hot_reg = zeros(1, 16);
      param_piston = zeros(1, 16);
      param_sphere = zeros(1, 4);

      % Declare external modules
      global mesh eos io bdry mesh_scaling;

      % Open file
      fid = fopen(filename, 'r');
      if fid == -1
        fprintf("ERROR: file %s doesn't exist.\n", filename);
        status = 1; % RC: returning status 1 might not be necessary
        return;
      end
      obj.filename = filename;

      % Read input
      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1); % RC: really doesn't seem to be used anywhere
      probname = fscanf(fid, '%s', 1);
      while any(ismember(inputtype(1), [' ', '#', '\', '$', '@', '&']))
        inputtype = fscanf(fid, '%s', 1);
        equalsigne = fscanf(fid, '%s', 1);
        probname = fscanf(fid, '%s', 1);
      end
      fprintf('%s : %s\n', inputtype, probname);

      if strcmp(probname, name_2shocks)
        inputtype = fscanf(fid, '%s', 1);
        equalsigne = fscanf(fid, '%s', 1);
        vel_move_in = fscanf(fid, '%f', 1);
        fprintf('%s : %e\n', inputtype, vel_move_in);
      end

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      dim = fscanf(fid, '%d', 1);
      fprintf('%s : %d\n', inputtype, dim);

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      ncell(1) = fscanf(fid, '%d', 1);
      fprintf('%s : %d\n', inputtype, ncell(1)); %this is x direction
      
      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      ncell(2) = fscanf(fid, '%d', 1);
      fprintf('%s : %d\n', inputtype, ncell(2)); % This is y

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      ncell(3) = fscanf(fid, '%d', 1);
      fprintf('%s : %d\n', inputtype, ncell(3)); % this is z

      %% mesh scaling
      % global mesh_scaling;
      ncell = round(2^mesh_scaling*ncell);
      if mod(ncell,2)~=0
        ncell = ncell + 1;
      end
      
      %%

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      nmat = fscanf(fid, '%d', 1);
      fprintf('%s : %d\n', inputtype, nmat);
      assert(nmat > 0 && nmat <= nmat_mx);

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      tmax_val = fscanf(fid, '%lf', 1);
      tmax = tmax_val;
      fprintf('%s : %e\n', inputtype, tmax_val);

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      ncycle_val = fscanf(fid, '%d', 1);
      ncycle = ncycle_val;
      fprintf('%s : %d\n', inputtype, ncycle_val);

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      dt_viz_val = fscanf(fid, '%lf', 1);
      dt_viz = dt_viz_val;
      fprintf('%s : %e\n', inputtype, dt_viz_val);

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      ncycle_viz_val = fscanf(fid, '%d', 1);
      ncycle_viz = ncycle_viz_val;
      fprintf('%s : %d\n', inputtype, ncycle_viz_val);

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      courant_val = fscanf(fid, '%lf', 1);
      courant = courant_val;
      fprintf('%s : %e\n', inputtype, courant_val);

      inputtype = fscanf(fid, '%s', 1);
      equalsigne = fscanf(fid, '%s', 1);
      dt_val = fscanf(fid, '%lf', 1);
      dt = dt_val;
      dt_initial = dt;
      fprintf('%s : %e\n', inputtype, dt_val);

      % Allocate materials and properties
      matids_ea_mat = 1:nmat;
      is_solid = zeros(1, nmat);   % logical array
      gamma_ea_mat = zeros(1, nmat);

      % Initialize based on problem name
      if strcmp(probname, name_accelpiston)
        disp("--- Piston problem -----------------------------");
        % Initialization for accelpiston
        for m = 1:nmat
          gamma_ea_mat(m) = 1.4;
        end
        nreg = 3;
        is_solid(1:nreg) = [0, 0, 1];
        reg2matids(1:nreg) = 1:nreg;
        rho_ea_reg(1:nreg) = [0.0, 1.0, 4.510];
        pres_ea_reg(1:nreg) = [0.0, 5.0e+09, 0.0];

        xl(1) = -50.0;
        xr(1) = 50.0;
        xl(2) = 0.0;
        xr(2) = 10.0;
        xl(3) = 0.0;
        xr(3) = 10.0;

        if dim == 2
          param_left(1:4) = [-50.0, 0.0, 0.0, 10.0];
        elseif dim == 3
          param_left(1:6) = [-50.0, 0.0, 0.0, 0.0, 10.0, 10.0];
        end
        reg_shape(2) = region_shape(shape_type.shape_rectangular, param_left);

        if dim == 2
          param_block(1:4) = [0.0, 0.0, 5.0, 10.0];
        elseif dim == 3
          param_block(1:6) = [0.0, 0.0, 0.0, 5.0, 10.0, 10.0];
        end
        reg_shape(3) = region_shape(shape_type.shape_rectangular, param_block);

        if dim == 2
          param_sphere(1:3) = [2.5, 5.0, 5.0];
        elseif dim == 3
          param_sphere(1:4) = [2.5, 5.0, 5.0, 5.0];
        end
        reg_shape(4) = region_shape(shape_type.shape_sphere, param_sphere);

        for reg = 1:nreg
          v_ea_reg(reg, 1:dim) = 0.0;
        end

      elseif strcmp(probname, name_shockformation)
        % Initialization for shockformation
        for m = 1:nmat
          gamma_ea_mat(m) = 1.4;
        end
        nreg = 2;
        reg2matids(1:nreg) = min(1:(nreg), nmat); %REVISIT
        rho_ea_reg(1:nreg) = [1.0e-06, 1.0e+02];
        pres_ea_reg(1:nreg) = [1.0e+04, 1.0e+04];

        xl(1) = 0.0;
        xr(1) = 50.0;
        ncell(1) = 500;
        xl(2:3) = 0.0;
        xr(2:3) = 1.0;
        ncell(2:3) = 10;

        if dim == 2
          param_left(1:4) = [0.0, 0.0, 0.2, 1.0];
        elseif dim == 3
          param_left(1:6) = [0.0, 0.0, 0.0, 0.2, 1.0, 1.0];
        end
        reg_shape(2) = region_shape(shape_type.shape_rectangular, param_left);

        for reg = 1:nreg
          v_ea_reg(reg, 1:dim) = 0.0;
        end
        v_ea_reg(2, 1) = 1.18e+05;

      elseif strcmp(probname, name_blowoff)
        % Initialization for blowoff
        for m = 1:nmat
          gamma_ea_mat(m) = 1.667;
        end
        nreg = 2;
        reg2matids(1:nreg) = min(1:(nreg), nmat); %REVISIT
        rho_ea_reg(1:nreg) = [0.0, 0.20];
        pres_ea_reg(1:nreg) = [0.0, 2.4e+10];

        xl(1) = -25.0;
        xr(1) = 25.0;
        ncell(1) = 500;
        xl(2:3) = 0.0;
        xr(2:3) = 1.0;
        ncell(2:3) = 10;

        if dim == 2
          param_left(1:4) = [-25.0, 0.0, 0.0, 1.0];
        elseif dim == 3
          param_left(1:6) = [-25.0, 0.0, 0.0, 0.0, 1.0, 1.0];
        end
        reg_shape(2) = region_shape(shape_type.shape_rectangular, param_left);

        for reg = 1:nreg
          v_ea_reg(reg, 1:dim) = 0.0;
        end

      elseif strcmp(probname, name_advection)
        % Initialization for advection
        gamma_ea_mat(1) = 1.4;
        if nmat > 1
          gamma_ea_mat(2) = 1.6;
        end

        nreg = 2;
        reg2matids(1) = 1;
        if nmat == 2 % REVISIT
          reg2matids(2) = 1;
        else
          for reg = 2:nreg % REVISIT
            m = reg;
            reg2matids(reg) = matids_ea_mat(m);
          end
        end
        rho_ea_reg(1:nreg) = [1.0, 10.0];
        pres_ea_reg(1:nreg) = [1.0, 1.0];

        if dim == 2
          xl(1:2) = [0.0, 0.0];
          xr(1:2) = [1.0, 0.2];
          ncell(1:2) = [100, 20];
        elseif dim == 3
          xl(1:3) = [0.0, 0.0, 0.0];
          xr(1:3) = [1.0, 0.2, 0.2];
          ncell(1:3) = [100, 20, 20];
        end

        if dim == 2
          param_left(1:4) = [0.0, 0.0, 0.1, 0.2];
        elseif dim == 3
          param_left(1:6) = [0.0, 0.0, 0.0, 0.1, 0.2, 0.2];
        end
        reg_shape(2) = region_shape(shape_type.shape_rectangular, param_left);

        for reg = 1:nreg
          v_ea_reg(reg, 1:dim) = 0.0;
          v_ea_reg(reg, 1) = 1.0;
        end

      elseif strcmp(probname, name_2shocks)
        % Initialization for 2shocks
        % gamma_ea_mat(1) = 1.4;
        % if nmat > 1
        %   gamma_ea_mat(2) = 1.6;
        % end
        %
        % nreg = 2;
        % reg2matids(1) = 1;
        % if nmat == 2 %REVISIT
        %   reg2matids(2) = 1;
        % else
        %   for reg = 2:nreg
        %     m = reg;
        %     reg2matids(reg) = matids_ea_mat(m);
        %   end
        % end
        % rho_ea_reg(1:2) = [1.0, 1.0];
        % pres_ea_reg(1:2) = [1.0, 1.0];
        %
        % if dim == 2
        %   xl(1:2) = [-0.5, 0.0];
        %   xr(1:2) = [0.5, 0.2];
        %   % ncell(1:2) = [100, 20];
        % elseif dim == 3
        %   xl(1:3) = [-0.5, 0.0, 0.0];
        %   xr(1:3) = [0.5, 0.2, 0.2];
        %   % ncell(1:3) = [100, 20, 20];
        % end
        %
        % if dim == 2
        %   param_left(1:4) = [-0.5, 0.0, 0.0, 0.2];
        % elseif dim == 3
        %   param_left(1:6) = [-0.5, 0.0, 0.0, 0.0, 0.2, 0.2];
        % end
        % reg_shape(2) = region_shape(shape_type.shape_rectangular, param_left);
        %
        % for reg = 1:nreg
        %   v_ea_reg(reg, 1:dim) = 0.0;
        % end
        % v_ea_reg(1, 1) = -vel_move_in;
        % v_ea_reg(2, 1) = vel_move_in;
        %%
        nx = ncell(1);ny = ncell(2);nz = ncell(3);
        % Assign gamma values
        for m = 1:nmat
          gamma_ea_mat(m) = 1.4;
        end
        gamma_ea_mat(1) = 1.4;
        gamma_ea_mat(2) = 1.6;

        % Number of regions
        nreg = 2;

        % Assign region-to-material mapping
        reg2matids = zeros(1, nreg); % Preallocate
        reg2matids(1) = 1;  % MATLAB is 1-based (C index 0)

        if nmat == 1
          reg2matids(2) = 1;  % Both regions use material 1
        else
          for reg = 2:nreg
            m = reg;
            reg2matids(reg) = matids_ea_mat(m);  % matid_ea_mat is assumed 1-based
          end
        end

        % Initial density and pressure per region
        rho_ea_reg = zeros(1, nreg);
        pres_ea_reg = zeros(1, nreg);

        rho_ea_reg(1) = 1.0;
        pres_ea_reg(1) = 1.0;

        rho_ea_reg(2) = 1.0;
        pres_ea_reg(2) = 1.0;

        % Set domain bounds and mesh resolution
        xl = zeros(1, 3);
        xr = zeros(1, 3);

        if dim == 2
          xl(1) = -0.5;
          xr(1) =  0.5;
          dx = (xr(1) - xl(1)) / nx;

          xl(2) = 0.0;
          xr(2) = xl(2) + dx * ny;

          ncell(1) = nx;
          ncell(2) = ny;
          ncell(3) = 1;
        elseif dim == 3
          xl(1) = -0.5;
          xr(1) =  0.5;
          dx = (xr(1) - xl(1)) / nx;

          xl(2) = 0.0;
          xr(2) = xl(2) + dx * ny;

          xl(3) = 0.0;
          xr(3) = xl(3) + dx * nz;
        end

        % Define region shape
        param_left = zeros(1, 6); % Supports both 2D and 3D

        if dim == 2
          param_left(1) = -0.5;  % x0
          param_left(2) = 0.0;   % y0
          param_left(3) = 0.0;   % x1
          param_left(4) = 0.2;   % y1
        elseif dim == 3
          param_left(1) = -0.5;  % x0
          param_left(2) = 0.0;   % y0
          param_left(3) = 0.0;   % z0
          param_left(4) = 0.0;   % x1
          param_left(5) = 0.2;   % y1
          param_left(6) = 0.2;   % z1
        end

        reg_shape = repmat(region_shape(shape_type.shape_universe, []), 1, nreg);
        reg_shape(2) = region_shape(shape_type.shape_rectangular, param_left);

        % Velocity field
        v_ea_reg = zeros(nreg, 3);
        v_ea_reg(1, 1) = -vel_move_in;
        v_ea_reg(2, 1) =  vel_move_in;

        %%

      elseif strcmp(probname, name_sod)
        % Initialization for sod
        disp("--- Sod shock problem -----------------------------");
        for m = 1:nmat
          gamma_ea_mat(m) = 1.4;
        end
        nreg = 2;
        reg2matids(1) = 1;
        reg2matids(2) = nmat;

        rho_ea_reg(1:2) = [0.125, 1.0];
        pres_ea_reg(1:2) = [0.1, 1.0];
        dx = 1.0/ncell(1);


        if dim == 2
          xl(1:2) = [0.0, 0.0];
          xr(1:2) = [1.0, xl(2) + dx*ncell(2)];

          param_left(1:4) = [0.0, 0.0, ...
            0.5*(xl(1)+xr(1)), xr(2)];
        elseif dim == 3
          xl(1:3) = [0.0, 0.0, 0.0];
          xr(1:3) = [1.0, xl(2) + dx*ncell(2), xl(3) + dx*ncell(3)];

          param_left(1:6) = [0.0, 0.0, 0.0, ...
            0.5*(xl(1)+xr(1)), xr(2), xr(3)];
        end
        reg_shape(2) = region_shape(shape_type.shape_rectangular, param_left);

        for reg = 1:nreg
          v_ea_reg(reg, 1:dim) = 0.0;
        end

      elseif strcmp(probname, name_mach60)
        % Initialization for mach60
        gamma_ea_mat(1) = 1.4;

        nreg = 2;
        reg2matids(1:2) = 1:nreg; %REVISIT, Will's code says up to nmat

        rho_ea_reg(1) = 1.0;
        pres_ea_reg(1) = 0.1;
        rho_ea_reg(2) = 3.99666;
        pres_ea_reg(2) = 449.975;

        if dim == 2
          xl(1:2) = [0.0, 0.0];
          xr(1:2) = [1.0, 0.2];
          ncell(1:2) = [100, 20];

          param_left(1:4) = [0.0, 0.0, 0.405, 0.2];
        elseif dim == 3
          xl(1:3) = [0.0, 0.0, 0.0];
          xr(1:3) = [1.0, 0.2, 0.2];
          ncell(1:3) = [100, 20, 20];

          param_left(1:6) = [0.0, 0.0, 0.0, 0.405, 0.2, 0.2];
        end
        reg_shape(2) = region_shape(shape_type.shape_rectangular, param_left);

        for reg = 1:nreg
          v_ea_reg(reg, 1:dim) = 0.0;
        end
        v_ea_reg(1, 1) = -18.3661;
      end

      % Calculate internal energy and set mesh
      for reg = 1:nreg
        m = reg2matids(reg);
        if ~is_solid(m)
          ei_ea_reg(reg) = pres_ea_reg(reg) / (gamma_ea_mat(m) - 1.0);
        else
          ei_ea_reg(reg) = eos.e_mie_gruneisen(rho_ea_reg(reg), pres_ea_reg(reg));
        end
      end

     %% Set mesh
      % mesh.set_mesh(dim, xl, xr, ncell, nb, nmat,matids_ea_mat);
      obj.btype_lower(1:dim) = bdry_type.bdry_transmitted;
      obj.btype_upper(1:dim) = bdry_type.bdry_transmitted;
      % 
      % % Set materials
      % vel_ea_reg = arrayfun(@(x) v_ea_reg(x, :), 1:nreg, 'UniformOutput', false);
      % mesh.set_mesh_mat(dim, obj.btype_lower, obj.btype_upper, ...
      %   nmat, matids_ea_mat, is_solid, gamma_ea_mat, ...
      %   nreg, reg2matids, reg_shape, rho_ea_reg, ...
      %   pres_ea_reg, ei_ea_reg, vel_ea_reg);

      % Update problem bounds
      obj.xl_prob(1:dim) = xl(1:dim);
      obj.xr_prob(1:dim) = xr(1:dim);
      obj.dims = dim;
      obj.nbdry = nb;
      obj.ncells(1:dim) = ncell(1:dim);
      obj.nummat = nmat;

      obj.matids_ea_mat =  matids_ea_mat;
      obj.is_solid      =  is_solid;
      obj.gamma_ea_mat  =  gamma_ea_mat;
      obj.tmax          =  tmax;
      obj.ncycle        =  ncycle;
      obj.dt_viz        =  dt_viz;
      obj.ncycle_viz    =  ncycle_viz;
      obj.courant       =  courant;
      obj.dt            =  dt;
      obj.problemname   =  probname;

      % Close the file
      fclose(fid);
      status = 0;
    end


    function hello_world(obj, str)
      fprintf("%s\n",str);
    end
  end
end
