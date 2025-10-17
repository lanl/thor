classdef TTAdvection < handle
  % c_advection   advection methodsmesh

  properties
    small = 1e-8;
    % Intercell edge variables
    dvol_for_edge = [];
    dmass_for_edge = [];
    dener_for_edge = [];

    % Conserved quantities (4D  array, size [nz, ny, nx, nmat])
    vol_for_cell  = [];
    mass_for_cell = [];
    ener_for_cell = [];

  end
  methods

    function vel_face = tt_get_face_velocity(~, mesh, dir)

      % TENSORIZING
      dim = mesh.dim_prob;  % Dimensionality of the mesh
      ncell = mesh.ncell_prob;
      nbdry = mesh.nbdry_prob;

      % Compute extended cell and node sizes
      ncell_ext = ncell + 2 * nbdry;
      nnode_ext = ncell_ext + 1;

      if dim == 3
        if dir == 1
          K = 1:ncell_ext(3);
          J = 1:ncell_ext(2);
          I = 1:nnode_ext(1);
          % tt
          % interpolation in j-dimension
          ops = {'intp','intp','none'};
          Cidx = {{K,K+1},{J,J+1},{I,0}};
          ops_factor = [0.5,0.5,1.0];

        elseif dir == 2
          K = 1:ncell_ext(3);
          J = 1:nnode_ext(2);
          I = 1:ncell_ext(1);
          ops = {'intp','none','intp'};
          Cidx = {{K,K+1},{J,0},{I,I+1}};
          ops_factor = [0.5,1.0,0.5];
        elseif dir == 3
          K = 1:nnode_ext(3);
          J = 1:ncell_ext(2);
          I = 1:ncell_ext(1);
          ops = {'none','intp','intp'};
          Cidx = {{K,0},{J,J+1},{I,I+1}};
          ops_factor = [1.0,0.5,0.5];
        end
        vel_face = ...
          FD_tt_ops(mesh.vav_for_3dnode(:,:,:,dir), ops,Cidx, ops_factor);
        %%
      elseif dim == 2
        error('This file is for 3D !!!')
      end
    end

    function courant_adv = advection(obj, ttmesh, mat, ncycle, dir, dt, ...
        btype_lower, btype_upper, is_solid, gamma_ea_mat)

      % global tol c_checking h5_checking;
      %
      % global tt_tol;

      % dim = ttmesh.dim_prob;
      % ncell = ttmesh.ncell_prob;
      % nbdry = ttmesh.nbdry_prob;

      %% Get face velocity (single output)
      vel_facett = obj.tt_get_face_velocity(ttmesh, dir);

      %% mapping
      courant_adv = obj.mapping3d(dir, dt, ...
        btype_lower, btype_upper, is_solid, ...
        gamma_ea_mat, vel_facett, ncycle);

    end
    %% %%%% mapping 3D %%%%%%%%%%%%%%%%
    function courant_adv = mapping3d(obj, dir, dt, ...
        btype_lower, btype_upper, ...
        is_solid, gamma_ea_mat, vel_face_3d_tt, ncycle)

      global tiny h5_checking tol;

      %variable for tt
      global ttmesh tt_tol ttbdry ttadv;

      % Initialize variables
      dim = 3;
      ncell = ttmesh.ncell_prob;
      dx = ttmesh.dx;
      nbdry = ttmesh.nbdry_prob;
      nmat = ttmesh.nmat_mesh;

      % Compute grid dimension variables
      vol_cell = prod(dx);
      ncell_ext = ncell + 2*nbdry;
      ncell_bdry = ncell + nbdry; % cell at the right boundary
      nnode_ext = ncell_ext + 1;
      nnode_bdry = nnode_ext - nbdry; % node at the right boundary
      sizes_edge = ncell_ext;
      sizes_edge(dir) = nnode_ext(dir);


      %% Compute average density
      % tt
      temp = ttmesh.rho_for_3dcell(nbdry+1:ncell_bdry(3),...
        nbdry+1:ncell_bdry(2), nbdry+1:ncell_bdry(1));
      rho_average = sum(temp)/prod(temp.n);


      %% TT
      % compute courant
      dist = vel_face_3d_tt * dt;
      fracten = dist ./ dx(dir);
      courant_adv = max(0.0,tt_max_abs(2.0*fracten));

      %%
      obj.tt_quantities_crossing_edge_3d(dir,fracten);

      %% Step 3: Compute material update with tensorized function
      % tt
      obj.tt_update_materials_after_advection(dir, nmat, vol_cell, sizes_edge, ncycle);

      %% Mapping velocity
      % tt
      tempfun = @(k,j,i) cal_mom_for_node_3D(k,j,i,nnode_ext,ttmesh,vol_cell);
      mom_for_nodett = amen_cross_zero([nnode_ext(3),nnode_ext(2),nnode_ext(1)],...
        @(x) cross_fun_nD(x,tempfun,dim),tt_tol,'verb',0);
      % reshape
      mom_for_nodett = tt_reshape(mom_for_nodett,[nnode_ext(3),...
        nnode_ext(2),nnode_ext(1),dim]);

      %% tt - compute mom_for_3dnode_newtt
      tempfun = @(k,j,i) cal_mom_for_node_new_3D(k, j, i,...
        ttmesh.vel_for_3dnode, mom_for_nodett, dir, dt, dx, obj.small, dim, nnode_ext);
      mom_for_3dnode_newtt = amen_cross_zero([nnode_ext(3),nnode_ext(2),...
        nnode_ext(1)],@(x) cross_fun_nD(x,tempfun,dim),tt_tol,'verb',0);
      %reshape
      mom_for_3dnode_newtt = tt_reshape(mom_for_3dnode_newtt,[nnode_ext(3),nnode_ext(2),...
        nnode_ext(1),dim],tt_tol);


      %% Update cell variables for rho, ei, and pres based on the volume fractions

      %% tt
      %vf
      tempfun = @(k,j,i,m) update_vf_3dmat(k,j,i,m,ttmesh.vf_3dmat,...
        ttadv,vol_cell,nbdry,ncell_bdry);
      ttmesh.vf_3dmat = amen_cross_zero([ncell_ext(3),ncell_ext(2),...
        ncell_ext(1),nmat],@(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);

      %rho
      tempfun = @(k,j,i,m) update_rho_3dmat(k,j,i,m,ttmesh.rho_3dmat,ttadv,...
        nbdry,ncell_bdry,tiny);
      ttmesh.rho_3dmat = amen_cross_zero([ncell_ext(3),ncell_ext(2),...
        ncell_ext(1),nmat],@(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);

      %ei
      tempfun = @(k,j,i,m) update_ei_3dmat(k,j,i,m,ttmesh.ei_3dmat, ...
        ttadv, nbdry,ncell_bdry,tiny);
      ttmesh.ei_3dmat = amen_cross_zero([ncell_ext(3),ncell_ext(2),...
        ncell_ext(1),nmat],@(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);
      %pres
      tempfun = @(k,j,i,m) update_pres_3dmat(k,j,i,m,ttmesh,...
        nbdry,ncell_bdry,is_solid,gamma_ea_mat);
      ttmesh.pres_3dmat = amen_cross_zero([ncell_ext(3),ncell_ext(2),...
        ncell_ext(1),nmat],@(x) cross_fun_nD(x,tempfun),tt_tol,'verb',0);

      %% Aggregate material properties into cell-wide properties
      %% tt
      temp = round(ttmesh.vf_3dmat.*ttmesh.rho_3dmat,tt_tol);
      G = core2cell(temp);
      G{3} = tensorprod(G{3},sum(G{4},2),3,1);
      ttmesh.rho_for_3dcell = cell2core(tt_tensor,G(1:3));

      temp = round(ttmesh.vf_3dmat.*ttmesh.ei_3dmat,tt_tol);
      G = core2cell(temp);
      G{3} = tensorprod(G{3},sum(G{4},2),3,1);
      ttmesh.ei_for_3dcell = cell2core(tt_tensor,G(1:3));

      temp = round(ttmesh.vf_3dmat.*ttmesh.pres_3dmat,tt_tol);
      G = core2cell(temp);
      G{3} = tensorprod(G{3},sum(G{4},2),3,1);
      ttmesh.pres_for_3dcell = cell2core(tt_tensor,G(1:3));

      %% Apply boundary conditions to update the boundary cells
      %tt
      ttbdry.tt_bdry_cell_3d(ttmesh,btype_lower, btype_upper);

      %% tt
      tempfun = @(k,j,i) update_velocity(k,j,i,nbdry,nnode_bdry,ttmesh,...
        mom_for_3dnode_newtt,rho_average,vol_cell,dim,obj.small);

      temptt = amen_cross_zero([nnode_ext(3),nnode_ext(2),nnode_ext(1)],...
        @(x) cross_fun_nD(x,tempfun,dim),tt_tol,'verb',0);
      ttmesh.vel_for_3dnode = reshape(temptt,...
        [nnode_ext(3),nnode_ext(2),nnode_ext(1),dim]);

      %% Apply boundary conditions to update the boundary cells
      %tt
      ttbdry.tt_bdry_node_3d(ttmesh, btype_lower, btype_upper);

    end
    %%
    %%
    function tt_quantities_crossing_edge_3d(obj, dir,fracten)
      global ttmesh tt_tol mat;

      ncell = ttmesh.ncell_prob;
      nbdry = ttmesh.nbdry_prob;
      dx = ttmesh.dx;
      xl_prob = ttmesh.xl_prob;
      nmat = ttmesh.nmat_mesh;

      vol_cell = prod(dx);
      ncell_ext = ncell + 2*nbdry;
      ncell_bdry = ncell + nbdry;
      nnode_ext = ncell_ext + 1;
      nnode_bdry = nnode_ext - nbdry;

      sizes_edge = ncell_ext;
      sizes_edge(dir) = nnode_ext(dir);

      % Set bounds for loops
      imax = ncell_bdry(1);
      jmax = ncell_bdry(2);
      kmax = ncell_bdry(3);
      %%
      if dir == 1, imax = nnode_bdry(1); end
      if dir == 2, jmax = nnode_bdry(2); end
      if dir == 3, kmax = nnode_bdry(3); end

      %%
      tempfun = @(k,j,i) cal_vme_3D(k, j, i, ttmesh, mat, obj, dir, fracten, nbdry,...
        kmax, jmax, imax, nmat, vol_cell,dx, xl_prob);
      % ytemp = tt_rand([sizes_edge(3),sizes_edge(2),sizes_edge(1)],3,1);
      Xtemp = amen_cross_zero([sizes_edge(3),sizes_edge(2),sizes_edge(1)], ...
        @(x) cross_fun_nD(x,tempfun,6), tt_tol,'verb',0);
      %reshape to split variables
      Xtemp = tt_reshape(Xtemp,[sizes_edge(3),sizes_edge(2),sizes_edge(1),nmat,3]);

      % matid_for_edge_tensortt = round(Xtemp(:,:,:,1), tt_tol);
      obj.dvol_for_edge = round(Xtemp(:,:,:,:,1), tt_tol);
      obj.dmass_for_edge = round(Xtemp(:,:,:,:,2), tt_tol);
      obj.dener_for_edge = round(Xtemp(:,:,:,:,3), tt_tol);
    end %ten_quantities_crossing_edge

    %%
    function tt_update_materials_after_advection(obj, dir, nmat, vol_cell, sizes_edge, ncycle)
      %TEN_UPDATE_MATERIALS_AFTER_ADVECTION Update material properties after tensor advection

      % global mesh mat h5_checking tol;

      global ttmesh tt_tol mat;

      ncell = ttmesh.ncell_prob;
      nbdry = ttmesh.nbdry_prob;

      ncell_ext = ncell + 2*nbdry;
      ncell_bdry = ncell_ext - nbdry - 1;



      %%
      obj.vol_for_cell = round(vol_cell * ttmesh.vf_3dmat,tt_tol);
      obj.mass_for_cell = round(obj.vol_for_cell .* ttmesh.rho_3dmat,tt_tol);
      obj.ener_for_cell = round(obj.vol_for_cell .* ttmesh.ei_3dmat,tt_tol);
      %
      % tt_check_obj(obj,ttadv,'Initialize volumed based on vf, rho and ei')
      %% Compute cell centric advect
      % tt
      obj = tt_cell_centric_update_materials_after_advection_3d(obj, dir);
      % tt_check_obj(obj,ttadv,'cell centric advection')

      %% tt
      tempfun = @(k,j,i) take_out_small_vof_3D(k,j,i,obj, mat, nbdry, ...
        ncell_bdry, nmat, vol_cell);

      Xtemp = amen_cross_zero([ncell_ext(3), ncell_ext(2),ncell_ext(1)],...
        @(x) cross_fun_nD(x,tempfun,6),tt_tol,'verb',0);

      Xtemp = tt_reshape(Xtemp,[ncell_ext(3),ncell_ext(2),ncell_ext(1),2,3]);

      obj.vol_for_cell = Xtemp(:,:,:,:,1);
      obj.mass_for_cell = Xtemp(:,:,:,:,2);
      obj.ener_for_cell = Xtemp(:,:,:,:,3);

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
            vf_local = full(mesh.vf_3dmat(k, j, i, 1:nmat));
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

  end
end
