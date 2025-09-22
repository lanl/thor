function ttmesh = tt_scale_2shock_mesh(ttmesh,mesh,scaling)
global tt_tol;
ttmesh.dim_prob = mesh.dim_prob;
ttmesh.nbdry_prob = mesh.nbdry_prob;
ttmesh.ncell_prob = mesh.ncell_prob.*(2^scaling);
ttmesh.xl_prob = mesh.xl_prob;
ttmesh.xr_prob = mesh.xr_prob;
ttmesh.dx = mesh.dx./(2^scaling);
ttmesh.nmat_mesh = mesh.nmat_mesh;
ttmesh.matids_mesh = mesh.matids_mesh;


ncell_prob= ttmesh.ncell_prob+4;
%%
% rho_for_3dcell = ones(ncell_prob(3),ncell_prob(2),ncell_prob(1))*mesh.rho_for_3dcell(1,1,1);
ttmesh.rho_for_3dcell = mesh.rho_for_3dcell(1,1,1)*tt_ones([ncell_prob(3),ncell_prob(2),ncell_prob(1)]);
% check_tt_error(rho_for_3dcell,ttmesh.rho_for_3dcell)
%%
% pres_for_3dcell = ones(ncell_prob(3),ncell_prob(2),ncell_prob(1))*mesh.pres_for_3dcell(1,1,1);
ttmesh.pres_for_3dcell = mesh.pres_for_3dcell(1,1,1)*...
  tt_ones([ncell_prob(3),ncell_prob(2),ncell_prob(1)]);
% check_tt_error(pres_for_3dcell,ttmesh.pres_for_3dcell)
%%
% ei_for_3dcell = ones(ncell_prob(3),ncell_prob(2),ncell_prob(1));
% ei_for_3dcell(:,:,1:ncell_prob(1)/2) = mesh.ei_for_3dcell(1,1,1);
% ei_for_3dcell(:,:,ncell_prob(1)/2+1:end) = mesh.ei_for_3dcell(1,1,end);

temp1 = tt_ones([ncell_prob(3),ncell_prob(2),ncell_prob(1)]);
G = core2cell(temp1);

G{3}(:,ncell_prob(1)/2+1:end,:) = 0;
G{3} = mesh.ei_for_3dcell(1,1,1)*G{3};
temp1 = cell2core(tt_tensor,G);

temp2 = tt_ones([ncell_prob(3),ncell_prob(2),ncell_prob(1)]);
G = core2cell(temp2);
G{3}(:,1:ncell_prob(1)/2,:) = 0;
G{3} = mesh.ei_for_3dcell(1,1,end)*G{3};
temp2 = cell2core(tt_tensor,G);
ttmesh.ei_for_3dcell = round(temp1 + temp2, tt_tol);

% check_tt_error(ei_for_3dcell,ttmesh.ei_for_3dcell)
%%
% vel_for_3dnode = zeros(ncell_prob(3)+1,ncell_prob(2)+1,ncell_prob(1)+1,3);
% vel_for_3dnode(:,:,1:ncell_prob(1)/2,1) = 1;
% vel_for_3dnode(:,:,ncell_prob(1)/2+2:end,1) = -1;

G = cell(1,4);
G{1} = ones(1,ncell_prob(3)+1,1);
G{2} = ones(1,ncell_prob(2)+1,1);
G{3} = zeros(1,ncell_prob(1)+1,1);
G{3}(:,1:ncell_prob(1)/2,:) = 1;
G{3}(:,ncell_prob(1)/2+2:end,:) = -1;
G{4} = [1,0,0];

ttmesh.vel_for_3dnode = cell2core(tt_tensor,G);


% check_tt_error(vel_for_3dnode,ttmesh.vel_for_3dnode)

%%
% vav_for_3dnode = zeros(ncell_prob(3)+1,ncell_prob(2)+1,ncell_prob(1)+1,3);
ttmesh.vav_for_3dnode = tt_zeros([ncell_prob(3)+1,ncell_prob(2)+1,ncell_prob(1)+1,3]);

%%
% cs_for_3dcell = ones(ncell_prob(3),ncell_prob(2),ncell_prob(1))*mesh.cs_for_3dcell(1,1,1);
% cs_for_3dcell(3:end,3:end,3:ncell_prob(1)/2) = mesh.cs_for_3dcell(3,3,3);

%% tt
tempfun = @(k,j,i) compute_cs_for_3dcell(k, j, i, mesh, ncell_prob);

ttmesh.cs_for_3dcell = amen_cross_zero([ncell_prob(3),ncell_prob(2),ncell_prob(1)], ...
  @(x) cross_fun_nD(x,tempfun,1),tt_tol,'verb',0);

% check_tt_error(cs_for_3dcell,ttmesh.cs_for_3dcell)
%%
% mesh.nmat_for_3dcell = ones(ncell_prob(3),ncell_prob(2),ncell_prob(1));
ttmesh.nmat_for_3dcell = tt_ones([ncell_prob(3),ncell_prob(2),ncell_prob(1)]);

%%
% matid_for_3dcell = ones(ncell_prob(3),ncell_prob(2),ncell_prob(1));
% matid_for_3dcell(:,:,1:ncell_prob(1)/2) = mesh.matid_for_3dcell(1,1,1);
% matid_for_3dcell(:,:,ncell_prob(1)/2+1:end) = mesh.matid_for_3dcell(1,1,end);

temp1 = tt_ones([ncell_prob(3),ncell_prob(2),ncell_prob(1)]);
G = core2cell(temp1);

G{3}(:,ncell_prob(1)/2+1:end,:) = 0;
G{3} = mesh.matid_for_3dcell(1,1,1)*G{3};
temp1 = cell2core(tt_tensor,G);

temp2 = tt_ones([ncell_prob(3),ncell_prob(2),ncell_prob(1)]);
G = core2cell(temp2);
G{3}(:,1:ncell_prob(1)/2,:) = 0;
G{3} = mesh.matid_for_3dcell(1,1,end)*G{3};
temp2 = cell2core(tt_tensor,G);


ttmesh.matid_for_3dcell = round(temp1 + temp2, tt_tol);

%%
% mesh.mixcell_for_3dcell = -1*ones(ncell_prob(3),ncell_prob(2),ncell_prob(1));

ttmesh.mixcell_for_3dcell = -1*tt_ones([ncell_prob(3),ncell_prob(2),ncell_prob(1)]);

%% compute mat variables
tempfun = @(k,j,i) compute_cell_properties(k, j, i, ttmesh);
temp = amen_cross_zero([ncell_prob(3),ncell_prob(2),ncell_prob(1)], ...
  @(x) cross_fun_nD(x,tempfun,8),tt_tol,'verb',0);

temp = tt_reshape(temp,[ncell_prob(3),ncell_prob(2),ncell_prob(1),8]);

G = core2cell(temp);

Gtemp = G; Gtemp{4} = Gtemp{4}(:,1:2,:);
ttmesh.rho_3dmat = round(cell2core(tt_tensor,Gtemp),tt_tol);
Gtemp = G; Gtemp{4} = Gtemp{4}(:,3:4,:);
ttmesh.ei_3dmat = round(cell2core(tt_tensor,Gtemp),tt_tol);
Gtemp = G; Gtemp{4} = Gtemp{4}(:,5:6,:);
ttmesh.pres_3dmat = round(cell2core(tt_tensor,Gtemp),tt_tol);
Gtemp = G; Gtemp{4} = Gtemp{4}(:,7:8,:);
ttmesh.vf_3dmat = round(cell2core(tt_tensor,Gtemp),tt_tol);
end

%% subroutines
function val = compute_cs_for_3dcell(k, j, i, mesh, ncell_prob)
% compute_cs_for_3dcell returns the value of cs_for_3dcell at (k,j,i)
% using the initialization logic:
% - default value from mesh.cs_for_3dcell(1,1,1)
% - override in region: k >= 3, j >= 3, 3 <= i <= ncell_prob(1)/2

% Default value
val = mesh.cs_for_3dcell(1,1,1);

% Check if (k,j,i) falls in the override region
if k >= 3 && j >= 3 && i >= 3 && i <= ncell_prob(1)/2
  val = mesh.cs_for_3dcell(3,3,3);
end
end
function vec8 = compute_cell_properties(k, j, i, ttmesh)
% compute_cell_properties returns [rho1, rho2, ei1, ei2, pres1, pres2, vf1, vf2]
% at position (k,j,i) based on mesh and smesh data

% Extract basic fields
rho  = ttmesh.rho_for_3dcell(k,j,i);
ei   = ttmesh.ei_for_3dcell(k,j,i);
pres = ttmesh.pres_for_3dcell(k,j,i);

% Get material ID at cell
m = int8(ttmesh.matid_for_3dcell(k,j,i)); % 1 or 2
m1 = (m == 1);
m2 = ~m1;

% Fill material-specific arrays
rho_1 = rho * m1;
rho_2 = rho * m2;

ei_1 = ei * m1;
ei_2 = ei * m2;

pres_1 = pres * m1;
pres_2 = pres * m2;

% VF logic: 1 if present, 0 if not
% vf_1 = 1 - (rho_2 > 0); % i.e., if mat 2 exists, vf_1 = 0
% vf_2 = double(rho_2 > 0);
if m1
  vf_1 = 1; vf_2 = 0;
else
  vf_1 = 0; vf_2 = 1;
end

% Output vector
vec8 = [rho_1, rho_2, ei_1, ei_2, pres_1, pres_2, vf_1, vf_2];
end
