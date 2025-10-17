function tt_solver(probname, mesh_scaling)

%% Determine input file based on problem name
if strcmp(probname, '2shocks')
  filename = 'input_2shocks_3d';
elseif strcmp(probname, 'sod')
  filename = 'input_sod_3d';
else
  error('Unsupported problem name: %s', probname);
end

%% Set global constants
global tiny small tol tol2;
tiny = 1.0e-30;
small = 1.0e-13;
tol = 1.0e-8;
tol2 = 1.0e-8;

%% Global classes/namespaces
global init geom mesh mat io bdry eos cupdate util vof2d vof3d xio vof3d_lib;
geom = c_geom();
mesh = c_mesh();
io = c_io();
eos = c_eos();
util = c_util();
xio = xdmfio();
cupdate = c_update();
vof2d = c_vof2d();
vof3d = c_vof3d();
vof3d_lib = c_lib_vof3d();

%% Specialized classes
global smesh sbdry update smat;
smesh = SMesh3D();
smat = SMat();
sbdry = SBdry();
update = SUpdate();

%% TT initialization
global ttmesh ttbdry ttupdate tt_tol ttadv;
tt_tol = 1e-6;
ttmesh = TTMesh3D();
ttupdate = TTUpdate();
ttbdry = TTBdry();
ttadv = TTAdvection();

%% Run Initialization
[init, err] = SInit2(filename);
if err ~= 0
  error('Initialization failed with error code %d', err);
end

%% Assign material and boundary
global bdry mat;
bdry = sbdry;
mat = smat;
mat.vfmin = 1e-8;

%% Apply mesh scaling based on problem
if strcmp(probname, '2shocks')
  data = load("Init_2shock.mat");
  ttmesh = tt_scale_2shock_mesh(ttmesh, data.mesh, mesh_scaling);
elseif strcmp(probname, 'sod')
  data = load("Init_sod.mat");
  ttmesh = tt_scale_sod_mesh(ttmesh, data.mesh, mesh_scaling);
end

%% Print mesh size
disp(" ");
cycletxt = sprintf('Mesh Size = %d x %d x %d', ttmesh.ncell_prob);
printBoxedText(cycletxt, 20);

%% Run the solver
ttcontrol = TTControl();
ttcontrol.run(init);
end