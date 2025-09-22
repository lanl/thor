classdef TTUpdate < handle
  %C_UPDATE Class for updating mesh properties
  %   This class contains methods to compute and update various physical
  %   quantities on a mesh, such as divergence of velocity, viscous forces,
  %   and energy, among others.


  properties
    % Force arrays
    force_for_2dnode  % 2D force at nodes, initially an empty 2D array
    force_for_3dnode  % 3D force at nodes, initially an empty 3D array

    % Viscous force (qvis) arrays
    qvis_for_2dcell   % 2D viscous force at cells, initially an empty 2D array
    qvis_for_3dcell   % 3D viscous force at cells, initially an empty 3D array

    % Divergence of velocity (divu) arrays
    divu_for_2dcell   % 2D divergence of velocity at cells, initially an empty 2D array
    divu_for_3dcell   % 3D divergence of velocity at cells, initially an empty 3D array
  end

  methods

    %% compute_divu
    %tt
    function tt_compute_divu(obj, mesh)
      % Compute the divergence of the velocity field in TT format
      % Inputs:
      %   - obj: instance of c_update
      %   - mesh: mesh in TT format

      % Shortcuts
      global tt_tol;
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      dx = mesh.dx(1:dim);
      ncell_ext = ncell + 2 * nbdry;  % Extended cell count including boundaries

      if dim == 2
        %dvx-tt
        % interpolation in j-1st dim, differentiation in i-2nd dim
        % vel_for_2dnode is node value. divu is cell value
        Gu = core2cell(mesh.vel_for_2dnode); % 2d has 2cores
        Gu{2} = tensorprod(Gu{2},Gu{3}(:,1,:),3,1);

        dvx = FD_tt_ops(Gu(1:2),{'intp','diff'},...
          {{1:(ncell_ext(2)),2:ncell_ext(2)+1},...
          {1:(ncell_ext(1)),2:ncell_ext(1)+1}},[0.5,1/dx(1)]);

        %dvy-tt
        Gu = core2cell(mesh.vel_for_2dnode); % 2d has 2cores
        Gu{2} = tensorprod(Gu{2},Gu{3}(:,2,:),3,1);
        dvy = FD_tt_ops(Gu(1:2),{'diff','intp'},...
          {{1:(ncell_ext(2)),2:ncell_ext(2)+1},...
          {1:(ncell_ext(1)),2:ncell_ext(1)+1}},[1/dx(2),0.5]);

        obj.divu_for_2dcell = round(dvx + dvy,tt_tol);

      elseif dim == 3
        % error('has not been tensorized !!!')
        % dvx-tt
        Gu = core2cell(mesh.vel_for_3dnode); %
        Gu{3} = tensorprod(Gu{3},Gu{4}(:,1,:),3,1);%get the first component

        dvx = FD_tt_ops(Gu(1:3),{'intp','intp','diff'},...
          {{1:(ncell_ext(3)),2:ncell_ext(3)+1},...
          {1:(ncell_ext(2)),2:ncell_ext(2)+1},...
          {1:(ncell_ext(1)),2:ncell_ext(1)+1}},[0.5,0.5,1/dx(1)]);

        % dvy-tt
        Gu = core2cell(mesh.vel_for_3dnode); %
        Gu{3} = tensorprod(Gu{3},Gu{4}(:,2,:),3,1);%get the second component

        dvy = FD_tt_ops(Gu(1:3),{'intp','diff','intp'},...
          {{1:(ncell_ext(3)),2:ncell_ext(3)+1},...
          {1:(ncell_ext(2)),2:ncell_ext(2)+1},...
          {1:(ncell_ext(1)),2:ncell_ext(1)+1}},[0.5,1/dx(1),0.5]);

        % dvy-tt
        Gu = core2cell(mesh.vel_for_3dnode); %
        Gu{3} = tensorprod(Gu{3},Gu{4}(:,3,:),3,1);%get the second component

        dvz = FD_tt_ops(Gu(1:3),{'diff','intp','intp'},...
          {{1:(ncell_ext(3)),2:ncell_ext(3)+1},...
          {1:(ncell_ext(2)),2:ncell_ext(2)+1},...
          {1:(ncell_ext(1)),2:ncell_ext(1)+1}},[1/dx(1),0.5,0.5]);

        obj.divu_for_3dcell = round(dvx + dvy + dvz, tt_tol);
      end
    end
    %% compute_qvis
    %tt
    function tt_compute_qvis(obj, mesh)
      % TT function
      global tt_tol;
      % Shortcuts
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      dx = mesh.dx(1);
      ncell_ext = ncell + 2 * nbdry;

      % Viscosity parameters
      c1 = 1.0;
      c2 = 0.1;

      if dim == 2
        tempfun = @(x1,x2,x3) qvis_2Dfun(obj,x1,x2,x3,dx,c1,c2);
        obj.qvis_for_2dcell = amen_cross_zero({obj.divu_for_2dcell,...
          mesh.rho_for_2dcell,mesh.cs_for_2dcell},@(x) cross_fun_nD(x,tempfun),...
          tt_tol,'verb',0);
      elseif dim == 3
        tempfun = @(x1,x2,x3) qvis_3Dfun(obj,x1,x2,x3,dx,c1,c2);
        obj.qvis_for_3dcell = amen_cross_zero({obj.divu_for_3dcell,...
          mesh.rho_for_3dcell,mesh.cs_for_3dcell},@(x) cross_fun_nD(x,tempfun),...
          tt_tol,'verb',0);
      end
    end
    %% compute_force
    %fg
    function tt_compute_force(obj, mesh, direction)
      % Compute the force for the mesh
      % Inputs:
      %   - obj: instance of c_update
      %   - mesh: instance of the mesh class with properties defining the problem
      %   - direction: direction of force calculation (0, 1, or 2)

      global tiny tt_tol ; % global variable defined in main


      % Extract required mesh properties
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      dx = mesh.dx(1:dim);

      % Initialize variables
      nnode_ext = ncell + 2 * nbdry + 1; % Extended nodes in each dimension
      vol_of_corner = prod(0.5 * dx);   % Volume of a corner
      area_of_node = prod(dx) / dx(direction); % Area of a node face perpendicular to the direction

      % Check dimensionality and initialize force arrays
      if dim == 2
        % Initialize the force_for_2dnode array if it is empty
        % obj.force_for_2dnode = zeros(nnode_ext(2), nnode_ext(1));


        % Calculate force for 2D nodes
        % for j = 2:nnode_ext(2) - 1
        %   j0 = j - 1;
        %   for i = 2:nnode_ext(1) - 1
        %     i0 = i - 1;
        % Compute mass for the node
        % mass = vol_of_corner * (...
        %   mesh.rho_for_2dcell(j0, i0) + mesh.rho_for_2dcell(j0, i) + ...
        %   mesh.rho_for_2dcell(j, i0) + mesh.rho_for_2dcell(j, i)) + tiny;

        %tt
        J =  2:(nnode_ext(2) - 1);
        I = 2:(nnode_ext(1) - 1);
        G = core2cell(mesh.rho_for_2dcell);
        mass = FD_tt_ops(G, {'intp','intp'},{{J,J-1},{I,I-1}},[vol_of_corner,1]);
        mass = round(mass + tiny,tt_tol);
        %%
        if direction == 1
          ops_list = {'intp','diff'};
        elseif direction == 2
          ops_list = {'diff','intp'};
        end
        
        %% compute force
        temp1 = FD_tt_ops(mesh.pres_for_2dcell, ops_list, {{J,J-1},{I,I-1}},[1,1]);
        temp2 = FD_tt_ops(obj.qvis_for_2dcell, ops_list, {{J,J-1},{I,I-1}},[1,1]);
        obj.force_for_2dnode = round(temp1 + temp2, tt_tol);

        
        %% Apply scaling
        % obj.force_for_2dnode(j, i) = obj.force_for_2dnode(j, i) * (area_of_corner / mass);
        massinv = amen_cross_zero({mass}, @(x) 1./x, tt_tol, 'verb', 0);
        obj.force_for_2dnode = (0.5 * area_of_node * obj.force_for_2dnode ) .* massinv;
        
        %% padd zeros to the boundary
        obj.force_for_2dnode = tt_pad_zeros(obj.force_for_2dnode);

      elseif dim == 3
        K = 2:(nnode_ext(3) - 1);
        J = 2:(nnode_ext(2) - 1);
        I = 2:(nnode_ext(1) - 1);

        G = core2cell(mesh.rho_for_3dcell);
        mass = FD_tt_ops(G, {'intp','intp','intp'},...
          {{K,K-1},{J,J-1},{I,I-1}},[vol_of_corner,1,1]);
        mass = round(mass + tiny,tt_tol);
        %%
        if direction == 1
          ops_list = {'intp','intp','diff'};
        elseif direction == 2
          ops_list = {'intp','diff','intp'};
        elseif direction == 3
          ops_list = {'diff','intp','intp'};
        end
        
        %% compute force
        temp1 = FD_tt_ops(mesh.pres_for_3dcell, ops_list,...
          {{K,K-1},{J,J-1},{I,I-1}},[1,1,1]);
        temp2 = FD_tt_ops(obj.qvis_for_3dcell, ops_list, ...
          {{K,K-1},{J,J-1},{I,I-1}},[1,1,1]);

        obj.force_for_3dnode = round(temp1 + temp2, tt_tol);
        %% Apply scaling
        % obj.force_for_2dnode(j, i) = obj.force_for_2dnode(j, i) * (area_of_corner / mass);
        massinv = amen_cross_zero({mass}, @(x) 1./x, tt_tol, 'verb', 0);
        obj.force_for_3dnode = (0.25 * area_of_node * obj.force_for_3dnode ) .* massinv;
        
        %% padd zeros to the boundary
        obj.force_for_3dnode = tt_pad_zeros(obj.force_for_3dnode);
      end
    end

    %% update_vel_comp
    function tt_update_vel_comp(obj, ttmesh, dt, dir)
      global tt_tol;
      % Some shortcuts
      dim = ttmesh.dim_prob;
      ncell = ttmesh.ncell_prob;
      nbdry = ttmesh.nbdry_prob;
      hdt = 0.5 * dt;  % Half of the time step
      nnode_ext = ncell + 2 * nbdry + 1;  % Extended number of nodes

      % Check dimensionality and update velocity components
      if dim == 2
        % Update velocity for 2D nodes
        %fg
        % for j = 2:nnode_ext(2)-1
        %   for i = 2:nnode_ext(1)-1
        %     hdv = hdt * obj.force_for_2dnode(j, i);
        %     mesh.vav_for_2dnode(j, i, dir) = mesh.vel_for_2dnode(j, i, dir) + hdv;
        %     mesh.vel_for_2dnode(j, i, dir) = mesh.vel_for_2dnode(j, i, dir) + 2 * hdv;
        %   end
        % end
        hdv = hdt*tt_get_inner(obj.force_for_2dnode); % inner part
        % update vav_for_2dnoe
        vav = ttmesh.vav_for_2dnode(:,:,dir);
        vavbdry = round(vav - tt_set_zero_boundaries(vav),tt_tol);
        vavinner = round(tt_get_inner(ttmesh.vel_for_2dnode(:,:,dir)) + hdv,tt_tol);
        newvav = round(vavbdry + tt_pad_zeros(vavinner),tt_tol);
        if dir==1
          ttmesh.vav_for_2dnode = round(tt_blockwise({newvav;ttmesh.vav_for_2dnode(:,:,2)}),tt_tol);
        elseif dir==2
          ttmesh.vav_for_2dnode = round(tt_blockwise({ttmesh.vav_for_2dnode(:,:,1);newvav}),tt_tol);
        end

        %update vel
        vel = round(ttmesh.vel_for_2dnode(:,:,dir) + tt_pad_zeros(2*hdv),tt_tol);
        if dir==1
          ttmesh.vel_for_2dnode = round(tt_blockwise({vel;ttmesh.vel_for_2dnode(:,:,2)}),tt_tol);
        elseif dir==2
          ttmesh.vel_for_2dnode = round(tt_blockwise({ttmesh.vel_for_2dnode(:,:,1);vel}),tt_tol);
        end

      elseif dim == 3
        hdv = hdt*tt_get_inner(obj.force_for_3dnode); % inner part
        % update vav_for_2dnoe
        vav = ttmesh.vav_for_3dnode(:,:,:,dir);
        vavbdry = round(vav - tt_set_zero_boundaries(vav),tt_tol);
        vavinner = round(tt_get_inner(ttmesh.vel_for_3dnode(:,:,:,dir)) + hdv,tt_tol);
        newvav = round(vavbdry + tt_pad_zeros(vavinner),tt_tol);
        if dir==1
          ttmesh.vav_for_3dnode = round(tt_blockwise({newvav;...
            ttmesh.vav_for_3dnode(:,:,:,2); ttmesh.vav_for_3dnode(:,:,:,3)}),tt_tol);
        elseif dir==2
          ttmesh.vav_for_3dnode = round(tt_blockwise({ttmesh.vav_for_3dnode(:,:,:,1);...
            newvav; ttmesh.vav_for_3dnode(:,:,:,3)}),tt_tol);
        elseif dir==3
          ttmesh.vav_for_3dnode = round(tt_blockwise({ttmesh.vav_for_3dnode(:,:,:,1);...
            ttmesh.vav_for_3dnode(:,:,:,2); newvav}),tt_tol);
        end

        %update vel
        vel = round(ttmesh.vel_for_3dnode(:,:,:,dir) + tt_pad_zeros(2*hdv),tt_tol);
        if dir==1
          ttmesh.vel_for_3dnode = round(tt_blockwise({vel;...
            ttmesh.vel_for_3dnode(:,:,:,2);ttmesh.vel_for_3dnode(:,:,:,3)}),tt_tol);
        elseif dir==2
          ttmesh.vel_for_3dnode = round(tt_blockwise({ttmesh.vel_for_3dnode(:,:,:,1);...
            vel;ttmesh.vel_for_3dnode(:,:,:,3)}),tt_tol);
        elseif dir==3
          ttmesh.vel_for_3dnode = round(tt_blockwise({ttmesh.vel_for_3dnode(:,:,:,1);...
            ttmesh.vel_for_3dnode(:,:,:,2); vel}),tt_tol);
        end

      end
    end

    function cournt = tt_update_density(obj, mesh, dt)

      global tt_tol;
      % Some shorthands
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      cournt = 0.0;
      ncell_ext = ncell + 2 * nbdry;
      nmat = mesh.nmat_mesh;
      % Check dimensionality and update density
      if dim == 2

        tempfun = @(j,i,m) update_rho_2dmat_fun(j,i,m,obj,mesh, dt, ncell_ext);

        mesh.rho_2dmat = amen_cross_zero([ncell_ext(2),ncell_ext(1),nmat],...
          @(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);
        
        % update rho_for_2dcell
        tempfun = @(j,i) update_rho_for_2dcell_fun(obj,j,i,mesh,ncell_ext);

        mesh.rho_for_2dcell = amen_cross_zero([ncell_ext(2),ncell_ext(1)],...
          @(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);
        %compute cournt
        cournt = max(cournt,tt_max(-dt*obj.divu_for_2dcell));

        %%%%%%%
      elseif dim == 3
        tempfun = @(k,j,i,m) update_rho_3dmat_fun(k,j,i,m,obj,mesh, dt, ncell_ext);

        mesh.rho_3dmat = amen_cross_zero([ncell_ext(3), ncell_ext(2),ncell_ext(1),nmat],...
          @(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);

        % update rho_for_2dcell
        tempfun = @(k,j,i) update_rho_for_3dcell_fun(obj,k,j,i,mesh,ncell_ext);

        mesh.rho_for_3dcell = amen_cross_zero([ncell_ext(3),ncell_ext(2),ncell_ext(1)],...
          @(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);
        %compute cournt
        cournt = max(cournt,tt_max(-dt*obj.divu_for_3dcell));
      end

    end
    %%  %%% ENERGY %%%%%%%
    function tt_update_energy(obj, mesh, dt)

      global tiny tt_tol;

      % Some shorthands
      dim = mesh.dim_prob;
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;
      ncell_ext = ncell + 2 * nbdry;
      nmat = mesh.nmat_mesh;

      % Check dimensionality and update energy
      if dim == 2

        tempfun = @(j,i,m) update_energy2D_fun(j,i,m,obj,mesh,dt,ncell_ext,tiny);

        mesh.ei_2dmat = amen_cross_zero([ncell_ext(2),ncell_ext(1),nmat],...
          @(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);

      elseif dim == 3
        tempfun = @(k,j,i,m) update_energy3D_fun(k,j,i,m,obj,mesh,dt,ncell_ext,tiny);

        mesh.ei_3dmat = amen_cross_zero([ncell_ext(3), ncell_ext(2),ncell_ext(1),nmat],...
          @(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);

      end

    end
    
    %% %% PRESSURE
    function tt_update_pressure(obj, ttmesh, is_solid_ea_mat, gamma_ea_mat)
      global tt_tol mesh;
      dim = ttmesh.dim_prob;
      nbdry = ttmesh.nbdry_prob;
      ncell_ext = ttmesh.ncell_prob + 2 * nbdry;
      nmat = ttmesh.nmat_mesh;
      if dim == 2
       % mat variable
       tempfun = @(j,i,m) update_presmat2D_fun(obj,j,i,m,ttmesh,gamma_ea_mat,...
         is_solid_ea_mat);
       Xtemp = amen_cross_zero([ncell_ext(2),ncell_ext(1),nmat],...
          @(x) cross_fun_nD(x,tempfun,2),tt_tol,'verb',0);
       %separate tt
       G = core2cell(Xtemp);
       G1 = G;
       G1{3} = G{3}(:,:,1);
       ttmesh.ei_2dmat = cell2core(tt_tensor,G1);
       G1{3} = G{3}(:,:,2);
       ttmesh.pres_2dmat = cell2core(tt_tensor,G1);
       
       %cell variable
       tempfun = @(j,i) update_prescell2D_fun(obj,j,i,mesh);
       Xtemp = amen_cross([ncell_ext(2),ncell_ext(1)],...
          @(x) cross_fun_nD(x,tempfun,2),tt_tol,'verb',0);
       %separate tt
       Xtemp = tt_reshape(Xtemp,[ncell_ext(2),ncell_ext(1),2]);
       ttmesh.ei_for_2dcell = round(Xtemp(:,:,1),tt_tol);
       ttmesh.pres_for_2dcell = round(Xtemp(:,:,2),tt_tol);
 
      elseif dim == 3
        % mat variable
        tempfun = @(k,j,i,m) update_presmat3D_fun(obj,k,j,i,m,ttmesh,gamma_ea_mat,...
          is_solid_ea_mat);
        Xtemp = amen_cross_zero([ncell_ext(3), ncell_ext(2),ncell_ext(1),nmat],...
          @(x) cross_fun_nD(x,tempfun,2),tt_tol,'verb',0);
        %separate tt
        G = core2cell(Xtemp);
        G1 = G;
        G1{dim+1} = G{dim+1}(:,:,1);
        ttmesh.ei_3dmat = cell2core(tt_tensor,G1);
        G1{dim+1} = G{dim+1}(:,:,2);
        ttmesh.pres_3dmat = cell2core(tt_tensor,G1);

        %cell variable
        tempfun = @(k,j,i) update_prescell3D_fun(obj,k,j,i,mesh);
        Xtemp = amen_cross_zero([ncell_ext(3), ncell_ext(2),ncell_ext(1)],...
          @(x) cross_fun_nD(x,tempfun,2),tt_tol,'verb',0);

        %separate tt
        Xtemp = tt_reshape(Xtemp,[ncell_ext(3),ncell_ext(2),ncell_ext(1),2]);
        ttmesh.ei_for_3dcell = round(Xtemp(:,:,:,1),tt_tol);
        ttmesh.pres_for_3dcell = round(Xtemp(:,:,:,2),tt_tol);
      end
    end
    %% %%%%%%%%%%%% 
    % Routines for TT
    function qvis_for_2dcell = qvis_2Dfun(obj,divu_for_2dcell,rho_for_2dcell,...
        cs_for_2dcell,dx,c1,c2)
      if divu_for_2dcell < 0.0
        dv = dx * divu_for_2dcell;
        qvis_for_2dcell = abs(dv) * rho_for_2dcell * ...
          (c1 * cs_for_2dcell - c2 * dv);
      else
        qvis_for_2dcell = 0;
      end
    end
    %%%%%%%
    function qvis_for_3dcell = qvis_3Dfun(obj,divu_for_3dcell,rho_for_3dcell,...
        cs_for_3dcell,dx,c1,c2)
      if divu_for_3dcell < 0.0
        dv = dx * divu_for_3dcell;
        qvis_for_3dcell = abs(dv) * rho_for_3dcell * ...
          (c1 * cs_for_3dcell - c2 * dv);
      else
        qvis_for_3dcell = 0;
      end
    end

    %% %%%%%
    function rho_2dmat = update_rho_2dmat_fun(j,i,m,obj,mesh,dt,ncell_ext)
      
      % % Update Courant number
      % courant = max(courant, -obj.divu_for_2dcell(j, i) * dt);
      if j > 1 && j < ncell_ext(2) && i > 1 && i < ncell_ext(1)
        % Compute factor for updating density
        factor = 1.0 / (1.0 + obj.divu_for_2dcell(j, i) * dt);

        % Update material densities
        rho_2dmat = mesh.rho_2dmat(j, i, m) * factor;
      else
        rho_2dmat = mesh.rho_2dmat(j, i, m);
      end
    end
    function rho_3dmat = update_rho_3dmat_fun(k,j,i,m,obj,mesh,dt,ncell_ext)

      % % Update Courant number
      % courant = max(courant, -obj.divu_for_2dcell(j, i) * dt);
      if k > 1 && k < ncell_ext(2) && j > 1 && j < ncell_ext(2) && i > 1 && i < ncell_ext(1)
        % Compute factor for updating density
        factor = 1.0 / (1.0 + obj.divu_for_3dcell(k,j, i) * dt);

        % Update material densities
        rho_3dmat = mesh.rho_3dmat(k, j, i, m) * factor;
      else
        rho_3dmat = mesh.rho_3dmat(k, j, i, m);
      end
    end
    %%
    function rho_for_2dcell = update_rho_for_2dcell_fun(~,j,i,mesh,ncell_ext)

      % % Update Courant number
      % courant = max(courant, -obj.divu_for_2dcell(j, i) * dt);
      if j > 1 && j < ncell_ext(2) && i > 1 && i < ncell_ext(1)
        rho_for_2dcell= sum(mesh.vf_2dmat(j, i, :) .* mesh.rho_2dmat(j, i, :));

      else
        rho_for_2dcell = mesh.rho_for_2dcell(j, i);
      end
    end

    function rho_for_3dcell = update_rho_for_3dcell_fun(~,k,j,i,mesh,ncell_ext)

      % % Update Courant number
      % courant = max(courant, -obj.divu_for_2dcell(j, i) * dt);
      if k > 1 && k < ncell_ext(3) && j > 1 && j < ncell_ext(2) && i > 1 && i < ncell_ext(1)
        rho_for_3dcell= sum(mesh.vf_3dmat(k, j, i, :) .* mesh.rho_3dmat(k,j, i, :));

      else
        rho_for_3dcell = mesh.rho_for_3dcell(k,j, i);
      end
    end
    
    
    %% %%%
    function ei_2dmat = update_energy2D_fun(j,i,m,obj,mesh,dt,ncell_ext,tiny)
      % Update specific internal energy for all materials
      ei_2dmat = mesh.ei_2dmat(j,i,m);
      if ei_2dmat <1e-13
        ei_2dmat = 0;
      end
      ei_2dmat = ei_2dmat / (mesh.rho_2dmat(j,i,m) + tiny); % Specific internal energy

      % Adjust specific internal energy for interior cells considering pressure and divergence
      if j > 1 && j < ncell_ext(2) && i > 1 && i < ncell_ext(1)
        work = (mesh.pres_2dmat(j,i,m) / (mesh.rho_for_2dcell(j,i) + tiny)) ...
          * obj.divu_for_2dcell(j,i) * dt;
        ei_2dmat = ei_2dmat - work; % Update specific internal energy
      end

    end
    %% %%%
    function ei_3dmat = update_energy3D_fun(k,j,i,m,obj,mesh,dt,ncell_ext,tiny)
      % Update specific internal energy for all materials
      ei_3dmat = mesh.ei_3dmat(k,j,i,m);
      if ei_3dmat <1e-13
        ei_3dmat = 0;
      end
      ei_3dmat = ei_3dmat / (mesh.rho_3dmat(k,j,i,m) + tiny); % Specific internal energy

      % Adjust specific internal energy for interior cells considering pressure and divergence
      if  k > 1 && k < ncell_ext(3) && ...
          j > 1 && j < ncell_ext(2) && ...
          i > 1 && i < ncell_ext(1)
        work = (mesh.pres_3dmat(k,j,i,m) / (mesh.rho_for_3dcell(k,j,i) + tiny)) ...
          * obj.divu_for_3dcell(k,j,i) * dt;
        ei_3dmat = ei_3dmat - work; % Update specific internal energy
      end

    end

    %% Update pressure
    function X = update_presmat2D_fun(obj,j,i,m,mesh,...
        gamma_ea_mat,is_solid_ea_mat)
     % Update specific internal energy density
     ei_2dmat = mesh.ei_2dmat(j, i, m) * mesh.rho_2dmat(j, i, m);

     if ~is_solid_ea_mat(m)
       % Calculate pressure for non-solid materials
       pres_2dmat = (gamma_ea_mat(m) - 1.0) * ei_2dmat;
     else
       % Use Mie-Gruneisen EOS for solid materials
       pres_2dmat = p_mie_gruneisen(mesh.rho_2dmat(j, i, m), ...
         ei_2dmat);
     end
     % output
     X(1) = ei_2dmat;
     X(2) = pres_2dmat;

    end

    function X = update_presmat3D_fun(obj,k,j,i,m,mesh,...
        gamma_ea_mat,is_solid_ea_mat)
      % Update specific internal energy density
      ei_3dmat = mesh.ei_3dmat(k,j,i,m) * mesh.rho_3dmat(k,j,i,m);

      if ~is_solid_ea_mat(m)
        % Calculate pressure for non-solid materials
        pres_3dmat = (gamma_ea_mat(m) - 1.0) * ei_3dmat;
      else
        % Use Mie-Gruneisen EOS for solid materials
        pres_3dmat = p_mie_gruneisen(mesh.rho_3dmat(k,j, i, m), ...
          ei_3dmat);
      end
      % output
      X(1) = ei_3dmat;
      X(2) = pres_3dmat;

    end
    %%
    function X = update_prescell2D_fun(obj,j,i,mesh)
      ei_for_2dcell = 0.0;
      pres_for_2dcell = 0.0;
      
      for m = 1:mesh.nmat_mesh
        % Aggregate energy and pressure contributions from all materials
        ei_for_2dcell = ei_for_2dcell + ...
          mesh.vf_2dmat(j, i, m) * mesh.ei_2dmat(j, i, m);

        pres_for_2dcell = pres_for_2dcell + ...
          mesh.vf_2dmat(j, i, m) * mesh.pres_2dmat(j, i, m);
      end
      X(1) = ei_for_2dcell;
      X(2) = pres_for_2dcell;
    end
    %%
    function X = update_prescell3D_fun(obj,k,j,i,mesh)
      ei_for_3dcell = 0.0;
      pres_for_3dcell = 0.0;
      
      for m = 1:mesh.nmat_mesh
        % Aggregate energy and pressure contributions from all materials
        ei_for_3dcell = ei_for_3dcell + ...
          mesh.vf_3dmat(k, j, i, m) * mesh.ei_3dmat(k, j, i, m);

        pres_for_3dcell = pres_for_3dcell + ...
          mesh.vf_3dmat(k, j, i, m) * mesh.pres_3dmat(k, j, i, m);
      end
      X(1) = ei_for_3dcell;
      X(2) = pres_for_3dcell;
    end
  end
end