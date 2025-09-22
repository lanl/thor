classdef c_control < handle
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
        function run(obj, init)

            % Global variables
            global tiny debug;
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

            % Calculate cell sizes and extended number of cells
            lsize = 1;
            for i = 1:dim
                dx(i) = (xr_prob(i) - xl_prob(i)) / ncell(i);
                ncell_ext(i) = ncell(i) + 2 * nbdry;
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
        
            % Apply boundary conditions
            bdry.cell_bdry_condition(mesh, btype_lower, btype_upper);
            
            % % Pass the data from the material
            % [nmixcell, nmixcell_int, nmixcell_mpoly, nmat_in_mixcell, ijk_in_mixcell, ...
            %     matid_in_mixcell, vf_in_mixcell, rho_in_mixcell, pres_in_mixcell, ei_in_mixcell] = mat_pass_mix_data(obj);
            nmixcell = mat.nmixcell;
            t = 0;
            ncycle = 0;
            ncycle_to_dump = ncycle_viz_freq;
            t_to_dump = dt_viz_freq;
        
            fileid = -1;
            %viz_dump(probname, fileid, dim, xl_prob, xr_prob, ncell, nbdry, t, ncycle, nmat, matids, ...
            %    nmat_for_2dcell, matid_for_2dcell, rho_for_2dcell, ei_for_2dcell, pres_for_2dcell, vel_for_2dnode, ...
            %    nmat_for_3dcell, matid_for_3dcell, rho_for_3dcell, ei_for_3dcell, pres_for_3dcell, vel_for_3dnode);
        
            mesh.mesh_sound_speed(nmat, matids, is_solid, gamma_ea_mat);
            mat.mat_sound_speed(mesh, nmat, matids, is_solid, gamma_ea_mat);
        
            courant_cs = mesh.courant_from_cs(dt);
        
            if dim == 2
                for j = 1:ncell_ext(2)
                    for i = 1:ncell_ext(1)
                        for dir = 1:dim
                            courant_adv = abs(mesh.vel_for_2dnode(j, i, dir)) * dt / dx(dir);
                            courant_cs = max(courant_cs, courant_adv);
                        end
                    end
                end
            elseif dim == 3
                for k = 1:ncell_ext(3)
                    for j = 1:ncell_ext(2)
                        for i = 1:ncell_ext(1)
                            for dir = 1:dim
                                courant_adv = abs(mesh.vel_for_3dnode(k, j, i, dir)) * dt / dx(dir);
                                courant_cs = max(courant_cs, courant_adv);
                            end
                        end
                    end
                end
            end
        
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
        
            adv = c_advection();
            while (t < tmax) && (ncycle < ncycle_max)
                cournmx = 0;

                % Lagrangian phase
                update.compute_divu(mesh);

                update.compute_qvis(mesh);
        
                for dir = 1:dim
                    update.compute_force(mesh, dir);
                    update.update_vel_comp(mesh, dt, dir);
                end
                
                % Working arrays
                nmixcell = mat.nmixcell;
                nm_tot = 0;
                for mx = 1:nmixcell
                    nm_tot = nm_tot + mat.nmat_in_mixcell(mx);
                end
        
                % Copy the densities and energies in mixed cells from the
                % pervious step
                rho_in_mixcell_old = cell(nmixcell, 1);
                es_in_mixcell_old = cell(nmixcell, 1);
                if nmixcell > 0
                    for mx = 1:nmixcell
                        rho_in_mixcell_old{mx} = mat.rho_in_mixcell{mx};
                        es_in_mixcell_old{mx} = mat.ei_in_mixcell{mx} ...
                                                  ./ (mat.rho_in_mixcell{mx} + tiny);
                    end
                end
        
                if dim == 2
                    mesh.rho_for_2dcell_old(:) = mesh.rho_for_2dcell(:);
                    mesh.es_for_2dcell_old(:) = mesh.ei_for_2dcell(:) ...
                                            ./ (mesh.rho_for_2dcell(:) + tiny);
                elseif dim == 3
                    mesh.rho_for_3dcell_old(:) = mesh.rho_for_3dcell(:);
                    mesh.es_for_3dcell_old(:) = mesh.ei_for_3dcell(:) ...
                                            ./ (mesh.rho_for_3dcell(:) + tiny);
                end
        
                % Courant factor update
                courant_div = update.update_density(mesh, mat, dt);
                cournmx = max(cournmx, courant_div);
        
                update.update_energy(mesh, mat, dt, rho_in_mixcell_old, es_in_mixcell_old);

                % Free temporary arrays
                clear rho_in_mixcell_old es_in_mixcell_old;
        
                % Apply boundary conditions
                bdry.cell_bdry_condition(mesh, btype_lower, btype_upper)
 

if (debug && ncycle == -1)
     xio.write_dump(mesh, mat, 37707, 1);
end
                % Advection phase
                for dir = 1:dim
                    courant_adv = adv.advection(mesh, mat, ncycle, pass+1, dt, ...
                        btype_lower, btype_upper, is_solid, gamma_ea_mat);                    
                    cournmx = max(cournmx, courant_adv);
                    pass = mod(pass + 1, dim);
                end
if (debug && ncycle == -1)
     xio.write_dump(mesh, mat, 37708, 1);
     error("That's all folks!");
end
       
                if fileid >= 0
                    close_file(fileid);
                end
        
                pass_start = mod(pass_start + 1, dim);
                pass = pass_start;
        
                t = t + dt;
                ncycle = ncycle + 1;
        
                % Determine the next dt
                mesh.mesh_sound_speed(nmat, matids, is_solid, gamma_ea_mat);
                mat.mat_sound_speed(mesh, nmat, matids, is_solid, gamma_ea_mat);
                courant_cs = mesh.courant_from_cs(dt);
                
                cournmx = max(cournmx, courant_cs);
        
                fprintf(' ncycle = %d,  t = %e, dt = %e, courant = %e\n', ncycle, t, dt, cournmx);
        
                ncycle_to_dump = ncycle_to_dump - 1;
                t_to_dump = t_to_dump - dt;
                
                if (ncycle_to_dump <= 0) || (t_to_dump <= 0) || (t >= tmax)
                    xio.write_dump(mesh, mat, ncycle, t);
        
                    ncycle_to_dump = ncycle_viz_freq;
                    t_to_dump = dt_viz_freq;
                end
        
                % Determine the next dt
                if cournmx > courant
                    dt = dt * (courant / cournmx);
                elseif cournmx < courant
                    ddt = (courant / cournmx - 1.0) * dt;
                    dt = dt + (0.1 * ddt);
                end
        
                dt = min(dt, tmax - t);
            end
        
            if fileid >= 0
                close_file(fileid);
            end
        
            % Free allocated memory
            if dim == 2
                clear rho_for_2dcell_old es_for_2dcell_old obj.cs_for_2dcell
            elseif dim == 3
                clear rho_for_3dcell_old es_for_3dcell_old obj.cs_for_3dcell
            end
        
        end

    end
end