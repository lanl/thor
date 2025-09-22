classdef TTControl < handle
  %INIT_STRUCT Summary of this class goes here
  %   Detailed explanation goes here

  properties
    % todo
  end
  methods
    % Original C declaration:
    % void control(const char *probname,
    %      int dim, int *ncell, int nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper,
    %      double *xl_prob, double *xr_prob,
    %      int nmat, int *matids, int *is_solid, double *gamma_ea_mat,
    %      double dt_initial, double tmax, int ncycle_max, double courant,
    %      int ncycle_viz_freq, double dt_viz_freq)
    function [ncycle,ncell_ext] = run(obj,init)

      % Global variables
      global ttbdry ttmesh ttupdate ttadv;
      global mat
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

      % Initialize variables
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

      %% Calculate cell sizes and extended number of cells
      lsize = 1;
      for i = 1:dim
        dx(i) = (xr_prob(i) - xl_prob(i)) / ncell(i);
        ncell_ext(i) = ncell(i) + 2 * nbdry;
        nnode_ext(i) = ncell(i) + 2 * nbdry + 1;
        lsize = lsize * ncell_ext(i);
      end

      %%
      % nmixcell = mat.nmixcell;
      t = 0;
      ncycle = 0;
      fileid = -1;

      %% mesh sound speed
      %tt
      ttmesh.tt_mesh_sound_speed(nmat, matids, is_solid, gamma_ea_mat);

      %% Compute courant_cs
      courant_cs = ttmesh.tt_courant_from_cs(dt);
      courant_cs = max(courant_cs, max(max(max(courant_adv))));
      if dim==3
        G = core2cell(ttmesh.vel_for_3dnode);%
        for dir = 1:dim
          G{dim+1}(:,dir,:) = G{dim+1}(:,dir,:)*dt/dx(dir);
        end
        temptt = cell2core(tt_tensor,G);
        maxtemptt = tt_max_abs(temptt);
        courant_cs = max(courant_cs,maxtemptt);
      end
      %%
      %fg
      if courant_cs > courant
        dt = courant * dt_initial / courant_cs;
        fprintf('WARNING: dt_initial too big, used dt = %e\n', dt);
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
      dt0 = dt;

      %% time stepping
      elapsed_time = 0;
      while (t < tmax) && (ncycle < ncycle_max)
        step_starttime = datetime;
        dt = dt0;
        cournmx = 0;
        %% %%%%%%%% Lagrangian phase
        %divu
        ttupdate.tt_compute_divu(ttmesh); %tt
        %qvis
        % update.compute_qvis(mesh);%fg
        ttupdate.tt_compute_qvis(ttmesh);%tt

        % force and vel_comp
        for dir = 1:dim
          ttupdate.tt_compute_force(ttmesh, dir);
          ttupdate.tt_update_vel_comp(ttmesh, dt, dir);
        end
        %% Courant factor update
        % update energy
        %tt
        ttupdate.tt_update_energy(ttmesh,dt);
        %% update density and courant factor
        courant_div = ttupdate.tt_update_density(ttmesh,dt);
        cournmx = max(cournmx, courant_div);

        %% Update pressure
        %tt
        ttupdate.tt_update_pressure(ttmesh, is_solid, gamma_ea_mat);

        %%
        % Apply boundary conditions
        ttbdry.cell_bdry_condition(ttmesh,btype_lower, btype_upper);

        %% %%%%%%%%%%%%%% Advection phase %%%%%%%%%%%%%%%%%%%%%%
        for dir = 1:dim
          %fg
          courant_adv = ttadv.advection(ttmesh, mat, ncycle, pass+1, dt, ...
            btype_lower, btype_upper, is_solid, gamma_ea_mat);
          cournmx = max(cournmx, courant_adv);
          %%
          pass = mod(pass + 1, dim);
        end

        if fileid >= 0
          close_file(fileid);
        end

        pass_start = mod(pass_start + 1, dim);
        pass = pass_start;

        t = t + dt;
        ncycle = ncycle + 1;

        %% Determine the next dt
        ttmesh.tt_mesh_sound_speed(nmat, matids, is_solid, gamma_ea_mat);
        courant_cs = ttmesh.tt_courant_from_cs(dt);
        cournmx = max(cournmx, courant_cs);
        % Determine the next dt
        if cournmx > courant
          dt = dt * (courant / cournmx);
        elseif cournmx < courant
          ddt = (courant / cournmx - 1.0) * dt;
          dt = dt + (0.1 * ddt);
        end

        dt = min(dt, tmax - t);
        step_endtime = datetime;
        elapsed_time = elapsed_time + seconds(step_endtime-step_starttime);
        fprintf(' ncycle = %d,  t = %e, dt = %e, courant = %e\n', ncycle, t, dt, cournmx);
        % fprintf('elapsed time per cycle = %.2f \n', seconds(step_endtime-step_starttime));
      end

      if fileid >= 0
        close_file(fileid);
      end
      %% write out variable to compare
      global mesh_scaling;
      savefilename = sprintf('./%s_mesh_%d_TT.mat',probname,mesh_scaling);
      rho_tt = ttmesh.rho_for_3dcell;
      ei_tt = ttmesh.ei_for_3dcell;
      pres_tt = ttmesh.pres_for_3dcell;
      save(savefilename,'rho_tt','ei_tt','pres_tt','elapsed_time','ncycle');
    end %run
  end %methods
end