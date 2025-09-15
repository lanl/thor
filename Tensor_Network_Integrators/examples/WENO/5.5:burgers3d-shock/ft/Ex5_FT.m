%% Info
% Example 5.5 : 3D Burgers' equation in full-tensor format with a shock
% Developer   : M. Engin Danis
% Notes       :
% It is recommended to run this script on terminal:
% (1) do: export TTWENOLIB=/path/of/WENO-paper
% (2) do: matlab -nodesktop -nojvm -nosplash -nodisplay -singleCompThread -r "try, run('Ex5_FT.m'); catch err, disp(getReport(err)); end;exit;"
% This script writes 1d results (x, u_ft) for each grid level under "results" directory, which can be used for postprocessing if needed.
clear; clc; close all;

addpath(genpath('../../../../matlab/WENO/src/ft/'))
addpath(genpath('../../../../matlab/WENO/src/misc/'))
addpath(genpath('../../../../matlab/WENO/utils/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))

%% set number of grid levels 
num_lvls = 7;
 
%% Input parameters to ftFDWenoSolver.m (Nx, Ny, Nz will be set later)
Tfinal    = 1.5/pi;                               % Final Time of the simulation   
dtOption  = 1;                                    % Time Stepping Option           
CFL       = 0.5;                                  % CFL Number for adaptive time-stepping
Neq       = 1;                                    % Number of equations 
N1d       = 5*(2.^linspace(1,num_lvls,num_lvls)); % Number of cells in each direction
Lx        = 1;                                    % Domain length in the x-direction 
Ly        = 1;                                    % Domain length in the y-direction 
Lz        = 1;                                    % Domain length in the z-direction 
funFlux   = @F_burgers3d;                         % Function handle for flux vector computation in each direction 
funEig    = @Eig_burgers3d;                       % Function handle for computing eigenvalues of flux jacobian    
funInit   = @Q0;                                  % Function handle for initial conditions                        
funBC     = @BCperiodic;                          % Function handle for boundary conditions                       
isSource  = false;                                % Flag to use user-defined "sourceTerm.m" to handle source terms    
gam       = [];                                   % Ratio of specific heats (not used here, so it is left empty)
isConsErr = false;                                % Flag to calculate the local conservation error

%% For spead measurements, make sure to use only a single thread
maxNumCompThreads(1); 
fprintf("Using %d threads\n",maxNumCompThreads);

%% Setup library paths
%restoredefaultpath;
%run(getenv("TTWENOLIB")+"/ft_setup.m");
%% Setup equation specific path for flux and eigenvalue routines 
%addpath(fdwenolib.ft+"/equations/Burgers3d");

%% Setup arrays for speed measurement reporting
CPUtime   = zeros(num_lvls,1);

%% Set the path where results will be written. This directory will be created if needed
resultPath = "results";
if exist(resultPath, 'dir')
    cmd = "rm -rf "+resultPath;
    system(cmd);
end
mkdir(resultPath);

%% Run simulations for all grid levels
for lvl=1:length(N1d)
    %
    Nx = N1d(lvl);
    Ny = 1; % change this to Nx for at least 2d simulations
    Nz = 1; % change this to Nx for at least 2d simulations
    %
    fprintf("%s",repelem("-",70));fprintf("\n");
    fprintf("Starting to solve grid level=%d with Nx=%d Ny=%d Nz=%d\n",lvl,Nx,Ny,Nz);
    fprintf("%s",repelem("-",70));fprintf("\n");
    %
    data = ftFDWenoSolver(Tfinal, dtOption, CFL, Neq, Nx, Ny, Nz, Lx, Ly, Lz,...
                          funFlux, funEig, funInit, funBC, isSource, gam, isConsErr);
    % Set time spent by the solver
    CPUtime(lvl) = data.CPUtime;
    % Write data to file
    Q  = data.Q; % numerical result
    %
    fileID = fopen(sprintf("%s/u-lvl=%d.txt",resultPath,lvl),'w');
    %
    for i=1:Nx
        %
        fprintf(fileID,'%.15e\t%.15e\t\n',data.X(i,1,1),Q(1,i,1,1));
        %
    end
    fclose(fileID);
    %
    save_errors(resultPath+"/time.txt",CPUtime(1:lvl));
    %
    fprintf("\n");
end
%
fprintf("\n ------ ------------\n");
fprintf(" level | CPUtime (s)\n");
fprintf(" ------ ------------\n");
fprintf(" %5d | %10.4e\n",[(1:length(CPUtime))' CPUtime]');
fprintf("\n");