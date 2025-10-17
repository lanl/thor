function fg_sparse_solver(probname, mesh_scaling)

%% Determine input file based on problem name
if strcmp(probname, '2shocks')
  filename = 'input_2shocks_3d';
elseif strcmp(probname, 'sod')
  filename = 'input_sod_3d';
else
  error('Unsupported problem name: %s', probname);
end

%% Set global constants
global tiny small tol;
global mesh_scaling;

tiny = 1.0e-30;
small = 1.0e-13;
tol = 1.0e-8;

%% Global classes/namespaces
global init geom mesh mat io bdry eos cupdate util vof2d vof3d xio vof3d_lib;
geom = c_geom();
mesh = c_mesh();
mat = c_mat();
io = c_io();
bdry = c_bdry();  % Tensorizing
eos = c_eos();
util = c_util();
xio = xdmfio();
cupdate = c_update();
vof2d = c_vof2d();
vof3d = c_vof3d();
vof3d_lib = c_lib_vof3d();

%% Specialized classes
global smesh sbdry update smat;
smesh = SMesh();
smat = SMat();
sbdry = SBdry();
update = SUpdate();

%% Run initialization
[init, err] = SInit(filename);

if err
  fprintf("Initialization failed with error code: %d\n", err);
  return;
end

%% Convert mesh to sparse mesh
convert_mesh_to_smesh_fn(mesh, smesh);

% Update global references
mesh = smesh;
bdry = sbdry;
mat = smat;
mat.vfmin = 1e-8;

%% Run the main sparse solver
control = SControl();
control.run(init);

end