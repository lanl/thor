function main(filename)

% Check if there are enough arguments; if not, complain and exit
if nargin ~= 1
    fprintf('usage: main(inputFileName)\n');
    return;
end

% Define global constants
global tiny small tol debug;
tiny = 1.0e-30;
small = 1.0e-12;
tol = 1.0e-10;
debug = true;

global chk;
chk = 0;

% Global classes/namespaces
global init geom mesh mat io bdry eos update util vof2d vof3d xio;
geom = c_geom();
mesh = c_mesh();
mat = c_mat();
io = c_io();
bdry = c_bdry();
eos = c_eos();
update = c_update();
util = c_util();
vof2d = c_vof2d();
vof3d = c_vof3d();
xio = xdmfio();

% initialization
[init, err] = c_init(filename);

% Check for errors
if err
    printf("Error: %d\n", err);
    return;
end

% Main driver function
control = c_control();
control.run(init);
