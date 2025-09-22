clear;close all;

addpath(genpath('../../../matlab/utils/tt-toolbox/'))

%
rng(1);
maxNumCompThreads(1);
%
isA              = true;  % flag to determine if diffusion operator is used
isF              = true;  % flag to determine if convection operator is used
tp               = 15;    % time where prediction ends
th               = 10;    % time where training ends
tint_order       = 4;     % time integration order for training and prediction
Cdt              = 1;     % Factor multiplying the training time step to determine the prediction time step
skip             = 1;     % number of actual timesteps between two training time steps
eps_qtt          = 1e-4;
eps_tt           = 1e-4;  % tt truncation
gamma            = 5e9;   % regularization paramter
g1               = 1e5;
g2               = g1;
u                = [];    % @(xx) 1;function handle to determine the time dependent source term
%
dt = 0.02;
t  = 0:dt:15;
%
[~,K]=min(abs(t-th));
% 
K  =  K - tint_order;
%
%restoredefaultpath;
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

filename  = opinflib.datadir+"/cylinder/qavg.bin";
filenamep = opinflib.datadir+"/cylinder/qp.bin";
%
fid = fopen(filename,"r");
N   = fread(fid,1,"int");
Neq = fread(fid,1,"int");
Nt  = fread(fid,1,"int");
Q   = fread(fid,"double");
fclose(fid);
%
Q    = reshape(Q,[N Neq Nt]);
Q    = calcPrimary(Q);
Qavg = Q;
%% Mesh
meshInfo = opinflib.datadir+"/cylinder/mesh_connectivity_info/";

% load mesh data
vert = load(sprintf('%s/nodes.txt',meshInfo));
conn = load(sprintf('%s/connectivity.txt',meshInfo));
%
% Setup the mesh
%
r = -1/3;
s = -1/3;
%
x = zeros(N,1);
y = zeros(N,1);
%
for i = 1:N
  %
  p1 = conn(i,1);
  p2 = conn(i,2);
  p3 = conn(i,3);
  %
  x(i) = -0.5*(r+s)*vert(p1,1) + 0.5*(r+1)*vert(p2,1) + 0.5*(s+1)*vert(p3,1);
  y(i) = -0.5*(r+s)*vert(p1,2) + 0.5*(r+1)*vert(p2,2) + 0.5*(s+1)*vert(p3,2);
  %
end
%
[~,idxp]=min(abs(t-tp));
%
types={'qtt','rom','ttrom'};
%
res = cell(length(types),1);
for idx=1:length(types)
  %
  type = types{idx};
  %
  addpath(opinflib.(type));
  %
  fprintf("Processing: %s-Opinf\n",type)
  %
  if(type=="qtt")
    %
    res{idx} = opinf(Q,u,t,K,isA,isF,tint_order,dt,Cdt*dt,tp,gamma,eps_qtt);
    %
    res{idx}.xp_ft = full(res{idx}.xp,res{idx}.xp.n');
    %
    Qout = reshape(res{idx}.xp_ft,[N Neq]);
    %
  elseif(type=="rom")
    %
    idx_ttrom = find(strcmp(types, "ttrom"));
    if (isempty(idx_ttrom) || idx_ttrom>idx)
      idx_ttrom = -1;
    end
    %
    if(idx_ttrom==-1)
      rn = 33;
    else
      rn = res{idx_ttrom}.n;
    end
    %
    Qout = reshape(Q,[N*Neq Nt]);
    %
    res{idx} = opinf(Qout,u,t,K,rn,isA,isF,tint_order,dt,Cdt*dt,tp,g1,g2);
    %
    Qout = reshape(res{idx}.xp,[N Neq]);
    %
  elseif(type=="ttrom")
    %
    Qout = reshape(Q,[N*Neq Nt]);
    %
    cfg = struct();
    cfg.filename   = filenamep; 
    cfg.offset     = 12;        % bytes to skip header (3×4‐byte ints)
    cfg.targetGB   = 5;         % desired chunk size in GB
    cfg.r0         = 1;
    %cfg.crosstype  = "greedy2";
    cfg.crosstype  = "cross2d";
    %
    res{idx} = opinf(cfg,u,t,K,isA,isF,tint_order,dt,Cdt*dt,tp,g1,g2,eps_tt);
    %
    Qout = reshape(res{idx}.xp,[N Neq]);
    %
  else
    error("rom or qtt or ttrom");
  end

  str = evalc('disp(res{idx})');  % capture display output
  fid = fopen(sprintf('out-%s.txt',type), 'w');
  fprintf(fid, '%s', str);
  fclose(fid);
  rmpath(opinflib.(type));
  %
  write2vtk(conn, vert, Qout, ["Density","U","V","Pressure"], "opinf_"+type+".vtk", idxp);
  %
end
%
write2vtk(conn, vert, Qavg(:,:,idxp), ["Density","U","V","Pressure"], "sim.vtk", idxp);
%
for idx=1:length(types)
    fprintf("%s: %f\n",types{idx},res{idx}.Tpod)
end
%
function Qp = calcPrimary(Q)
  %
  Qp = zeros(size(Q));
  %
  rho  = Q(:,1,:);
  rhoU = Q(:,2,:);
  rhoV = Q(:,3,:);
  rhoE = Q(:,4,:);
  %
  U    = rhoU./rho;
  V    = rhoV./rho;
  Vel2 = U.^2 + V.^2;
  rhoe = rhoE - 0.5*rho.*Vel2;
  P    = 0.4*rhoe;
  %
  Qp(:,1,:) = rho;
  Qp(:,2,:) = U;
  Qp(:,3,:) = V;
  Qp(:,4,:) = P;
  %
end