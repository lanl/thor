classdef SControl < handle
  %INIT_STRUCT Summary of this class goes here
  %   Detailed explanation goes here

  properties
    % empty
  end
  methods
    function [ncycle,ncell_ext] = run(~, init)

      % Global variables
      global mesh bdry mat update xio;
      % Copy arguments
      probname = init.problemname;
      dim = init.dims;
      ncell = init.ncells;
      nbdry = init.nbdry;
      btype_lower = init.btype_lower;
      btype_upper = init.btype_upper;
      xl_prob = init.xl_prob;
      xr_prob = init.xr_prob;
      nmat = init.nummat;
      matids = init.matids_ea_mat;
      is_solid = init.is_solid;
      gamma_ea_mat = init.gamma_ea_mat;
      dt_initial = init.dt;
      tmax = init.tmax;
      ncycle_max = init.ncycle;
      courant = init.courant;
      ncycle_viz_freq = init.ncycle_viz;
      dt_viz_freq = init.dt_viz;

      %% Initialize variables
      filename = '';
      ncycle = 0;
      ncycle_to_dump = 0;
      t = 0;
      dt = dt_initial;
      ddt = 0;
      cournmx = 0;
      courant_cs = 0;
      courant_div = 0;
      courant_adv = 0;
      qvis_coef = 0.1;
      dx = zeros(1, 3);
      ncell_ext = zeros(1, 3);

      % Calculate cell sizes and extended number of cells
      lsize = 1;
      for i = 1:dim
        dx(i) = (xr_prob(i) - xl_prob(i)) / ncell(i);
        ncell_ext(i) = ncell(i) + 2 * nbdry;
        nnode_ext(i) = ncell(i) + 2 * nbdry + 1;

        lsize = lsize * ncell_ext(i);
      end

      if dim == 2
        % Allocate memory for 2D arrays
        mesh.cs_for_2dcell = zeros(ncell_ext(2), ncell_ext(1));
        mesh.rho_for_2dcell_old = zeros(ncell_ext(2), ncell_ext(1));
        mesh.es_for_2dcell_old = zeros(ncell_ext(2), ncell_ext(1));
      elseif dim == 3
        % Allocate memory for 3D arrays
        mesh.cs_for_3dcell = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));
        mesh.rho_for_3dcell_old = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));
        mesh.es_for_3dcell_old = zeros(ncell_ext(3), ncell_ext(2), ncell_ext(1));
      end

      nmixcell = mat.nmixcell;
      t = 0;
      ncycle = 0;
      ncycle_to_dump = ncycle_viz_freq;
      t_to_dump = dt_viz_freq;

      fileid = -1;
      mesh.mesh_sound_speed(nmat, matids, is_solid, gamma_ea_mat);

      courant_cs = mesh.courant_from_cs(dt);
      courant_adv = mesh.courant_from_vel(dt);
      cournmx = max(cournmx,courant_adv);


      if courant_cs > courant
        dt = courant * dt_initial / courant_cs;
        warning('dt_initial too big, used dt = %e\n', dt);
      else
        dt = dt_initial;
      end

      if dim == 2
        vartype = 'mio_face';
      else
        vartype = 'mio_zone';
      end

      pass_start = 0;
      pass = pass_start;

      adv = SAdvection();      dt0 = dt;

      %% time stepping
      elapsed_time = 0;
      while (t < tmax) && (ncycle < ncycle_max)
        step_starttime = datetime;
        dt = dt0;
        cournmx = 0;
        %% Lagrangian phase
        update.compute_divu(mesh);
        update.compute_qvis(mesh);

        for dir = 1:dim
          update.compute_force(mesh, dir);
          update.update_vel_comp(mesh, dt, dir);
        end

        %% Courant factor update
        % update energy
        update.update_energy(mesh, dt);

        %update density and courant factor
        courant_div = update.update_density(mesh, dt);
        cournmx = max(cournmx, courant_div);

        %% Update pressure
        update.update_pressure(mesh, is_solid, gamma_ea_mat);

        %%
        % Apply boundary conditions
        bdry.cell_bdry_condition(mesh, btype_lower, btype_upper)

        %% Advection phase


        for dir = 1:dim
          courant_adv = adv.advection(mesh, mat, ncycle, pass+1, dt, ...
            btype_lower, btype_upper, is_solid, gamma_ea_mat);
          cournmx = max(cournmx, courant_adv);
          pass = mod(pass + 1, dim);
        end

        if fileid >= 0
          close_file(fileid);
        end

        pass_start = mod(pass_start + 1, dim);
        pass = pass_start;

        % dt = round(dt,4);
        t = t + dt;
        ncycle = ncycle + 1;

        % Determine the next dt
        mesh.mesh_sound_speed(nmat, matids, is_solid, gamma_ea_mat);
        mat.mat_sound_speed(mesh, nmat, matids, is_solid, gamma_ea_mat);
        courant_cs = mesh.courant_from_cs(dt);

        cournmx = max(cournmx, courant_cs);
        courant_adv = mesh.courant_from_vel(dt);
        cournmx = max(cournmx,courant_adv);

        fprintf(' ncycle = %d,  t = %e, dt = %e, courant = %e\n', ncycle, t, dt, cournmx);

        %%
        ncycle_to_dump = ncycle_to_dump - 1;
        t_to_dump = t_to_dump - dt;

        % if (ncycle_to_dump <= 0) || (t_to_dump <= 0) || (t >= tmax)
        %   xio.write_dump(mesh, mat, ncycle, t);
        % 
        %   ncycle_to_dump = ncycle_viz_freq;
        %   t_to_dump = dt_viz_freq;
        % end

        % Determine the next dt
        if cournmx > courant
          dt = dt * (courant / cournmx);
        elseif cournmx < courant
          ddt = (courant / cournmx - 1.0) * dt;
          dt = dt + (0.1 * ddt);
        end

        dt = min(dt0, tmax - t);
        step_endtime = datetime;
        elapsed_time = elapsed_time + seconds(step_endtime-step_starttime);
      end % time stepping loop

      %%
      if fileid >= 0
        close_file(fileid);
      end

      global mesh_scaling;
      savefilename = sprintf('./%s_mesh_%d_FG.mat',probname,mesh_scaling);
      rho = mesh.rho_for_3dcell;
      ei = mesh.ei_for_3dcell;
      pres= mesh.pres_for_3dcell;
      save(savefilename,'rho','ei','pres','elapsed_time','ncycle');

    end

  end
end
