function convert_mesh_to_smesh_fn(mesh,smesh)
if mesh.dim_prob == 2
    nmat = numel(mesh.matids_mesh);
    [ny,nx] = size(mesh.rho_for_2dcell);
    properties = {
      'dim_prob', ...
      'nbdry_prob', ...
      'ncell_prob', ...
      'xl_prob', ...
      'xr_prob', ...
      'dx', ...
      'nmat_mesh', ...
      'matids_mesh', ...
      'rho_for_2dcell',...
      'ei_for_2dcell',...
      'pres_for_2dcell',...
      'vel_for_2dnode', ...
      'vav_for_2dnode', ...
      'nmat_for_2dcell', ...
      'matid_for_2dcell',...
      'mixcell_for_2dcell',...
      'cs_for_3dcell'} ;
    
    
    % Copy each property if it exists in both objects
    for k = 1:length(properties)
      propertyName = properties{k};
      if isprop(mesh, propertyName) && isprop(smesh, propertyName)
        smesh.(propertyName) = mesh.(propertyName);
      else
        warning('Property "%s" does not exist in one of the objects.', propertyName);
      end
    end
    
    % properties = {'rho_for_2dcell','ei_for_2dcell','pres_for_2dcell'};
    smesh.rho_2dmat = zeros(ny,nx,nmat);
    smesh.ei_2dmat = zeros(ny,nx,nmat);
    smesh.pres_2dmat = zeros(ny,nx,nmat);
    smesh.vf_2dmat = zeros(ny,nx,nmat);
    for j = 1:ny
      for i = 1:nx
          smesh.rho_2dmat(j,i,smesh.matid_for_2dcell(j,i)) = mesh.rho_for_2dcell(j,i);
          smesh.ei_2dmat(j,i,smesh.matid_for_2dcell(j,i)) = mesh.ei_for_2dcell(j,i);
          smesh.pres_2dmat(j,i,smesh.matid_for_2dcell(j,i)) = mesh.pres_for_2dcell(j,i);
      end
    end
    smesh.vf_2dmat = double((smesh.rho_2dmat>0));
    smesh.vf_2dmat(:,:,1) = 1 - sum(smesh.vf_2dmat(:,:,2:nmat),3);

else % dim_prob == 3

    nmat = numel(mesh.matids_mesh);
    [nz,ny,nx] = size(mesh.rho_for_3dcell);
    properties = {
      'dim_prob', ...
      'nbdry_prob', ...
      'ncell_prob', ...
      'xl_prob', ...
      'xr_prob', ...
      'dx', ...
      'nmat_mesh', ...
      'matids_mesh', ...
      'rho_for_3dcell',...
      'ei_for_3dcell',...
      'pres_for_3dcell',...
      'vel_for_3dnode', ...
      'vav_for_3dnode', ...
      'nmat_for_3dcell', ...
      'matid_for_3dcell',...
      'mixcell_for_3dcell',...
      'cs_for_3dcell'} ;
    
    
    % Copy each property if it exists in both objects
    for k = 1:length(properties)
      propertyName = properties{k};
      if isprop(mesh, propertyName) && isprop(smesh, propertyName)
        smesh.(propertyName) = mesh.(propertyName);
      else
        warning('Property "%s" does not exist in one of the objects.', propertyName);
      end
    end
    
    % 3-D arrays (one value per cell)
    rho  = mesh.rho_for_3dcell;   % [nz ny nx]
    ei   = mesh.ei_for_3dcell;    % [nz ny nx]
    pres = mesh.pres_for_3dcell;  % [nz ny nx]


    % properties = {'rho_for_2dcell','ei_for_2dcell','pres_for_2dcell'};

    %  smesh.rho_3dmat = [];
    % smesh.ei_3dmat = [];
    % smesh.pres_3dmat = [];
    % smesh.vf_3dmat = [];


    smesh.rho_3dmat = zeros(nz,ny,nx,nmat);
    smesh.ei_3dmat = zeros(nz,ny,nx,nmat);
    smesh.pres_3dmat = zeros(nz,ny,nx,nmat);
    smesh.vf_3dmat = zeros(nz,ny,nx,nmat);
    % for k = 1:nz
    %   for j = 1:ny
    %     for i = 1:nx
    %       k
    %       m = smesh.matid_for_3dcell(k,j,i);
    %       smesh.rho_3dmat(k,j,i,m) = mesh.rho_for_3dcell(k,j,i);
    %       smesh.ei_3dmat(k,j,i,m) = mesh.ei_for_3dcell(k,j,i);
    %       smesh.pres_3dmat(k,j,i,m) = mesh.pres_for_3dcell(k,j,i);
    %     end
    %   end
    % end

    % Logical masks
    matid = smesh.matid_for_3dcell;
    m1 = (matid == 1);
    m2 = ~m1;              % since only two materials (1 or 2)

    % Fill each material “slot”
    smesh.rho_3dmat(:,:,:,1)  = rho  .* m1;
    smesh.rho_3dmat(:,:,:,2)  = rho  .* m2;

    smesh.ei_3dmat(:,:,:,1)   = ei   .* m1;
    smesh.ei_3dmat(:,:,:,2)   = ei   .* m2;

    smesh.pres_3dmat(:,:,:,1) = pres .* m1;
    smesh.pres_3dmat(:,:,:,2) = pres .* m2;

    smesh.vf_3dmat(:,:,:,1) = double((smesh.matid_for_3dcell==1));
    smesh.vf_3dmat(:,:,:,2) = double((smesh.matid_for_3dcell==2));

    
    % smesh.vf_3dmat(:,:,:,1) = 1 - sum(smesh.vf_3dmat(:,:,:,2:nmat),4);

end % if dim_prob == 3    

end