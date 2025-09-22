clear;close all;

addpath(genpath('../../../matlab/utils/tt-toolbox/'))

%
maxNumCompThreads(1); 
opinf_types = {"rom","ttrom","ft","tt","qtt"};
%
isA              = true;  % flag to determine if diffusion operator is used
isF              = false; % flag to determine if convection operator is used
tp               = 0.9;   % time where prediction ends
th               = 0.25;  % time where training ends
tint_order       = 4;     % time integration order for training and prediction
Cdt              = 1;     % Factor multiplying the training time step to determine the prediction time step
skip             = 1;     % number of actual timesteps between two training time steps
eps_tt           = 1e-12; % tt truncation
gamma            = 0;     % regularization paramter
g_rom            = 1e-6;  % regularization paramter for rom
u                = [];    % @(xx) 1;function handle to determine the time dependent source term
retained_energy  = 1;
%
Lx = 1;
%
%restoredefaultpath;
%run(getenv("TTTOOLBOXHOME")+"/setup.m");
%run(getenv("TTTOOLBOXHOME")+"/setup.m");
%run(getenv("OPINFHOME")+"/setup.m");
%

%opinflib.basedir = "../../";
%
opinflib.datadir = "../../../matlab/data/thor-data-main";
%
opinflib.srcdir = "../../../matlab/tt-opinf/src";
%
opinflib.tt    = opinflib.srcdir + "/tt";
opinflib.qtt   = opinflib.srcdir + "/qtt";
opinflib.ft    = opinflib.srcdir + "/ft";
opinflib.rom   = opinflib.srcdir + "/rom";
opinflib.ttrom = opinflib.srcdir + "/ttrom";
%
addpath(opinflib.srcdir);

filename = opinflib.datadir+"/heat";
%
fun_grid      = @(Nx,Ny,Nz) fun_grid1d(Nx,Ny,Nz,Lx,Lx,Lx);
func_grid_idx = @(Nx,Ny,Nz) func_grid1d_idx(Nx,Ny,Nz);
%
res = postprocess(opinflib, filename, opinf_types, isA, isF, tp, th, tint_order, Cdt, skip, eps_tt, ...
                  gamma, g_rom, u, retained_energy, fun_grid, func_grid_idx);
%
function grid1d = fun_grid1d(Nx,Ny,Nz,Lx,Ly,Lz)
  %
  dx = Lx/Nx;
  %
  gridx = ((dx:dx:Lx)-dx/2)';
  %
  grid1d = sqrt(2*gridx.^2);
  %
end
%
function grid1didx = func_grid1d_idx(Nx,Ny,Nz)
  %
  grid1didx = 1:Nx+1:Nx^2;
  %
end