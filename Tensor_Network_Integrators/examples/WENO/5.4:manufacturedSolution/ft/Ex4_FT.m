%% Info
% Example 5.4 : Manufactured Solution in full-tensor format
% Developer   : M. Engin Danis
% Notes       :
% It is recommended to run this script on terminal:
% (1) do: export TTWENOLIB=/path/of/WENO-paper
% (2) do: matlab -nodesktop -nojvm -nosplash -nodisplay -singleCompThread -r "try, run('Ex4_FT.m'); catch err, disp(getReport(err)); end;exit;"
clear; clc; close all;

addpath(genpath('../../../../matlab/WENO/src/ft/'))
addpath(genpath('../../../../matlab/WENO/src/misc/'))
addpath(genpath('../../../../matlab/WENO/utils/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))

%% set number of grid levels 
num_lvls = 4;
 
%% Input parameters to ftFDWenoSolver.m (Nx, Ny, Nz will be set later)
Tfinal    = 0.1;                                  % Final Time of the simulation   
dtOption  = 2;                                    % Time Stepping Option           
CFL       = 0.5;                                  % CFL Number for adaptive time-stepping
Neq       = 5;                                    % Number of equations 
N1d       = 5*(2.^linspace(1,num_lvls,num_lvls)); % Number of cells in each direction
Lx        = 1;                                    % Domain length in the x-direction 
Ly        = 1;                                    % Domain length in the y-direction 
Lz        = 1;                                    % Domain length in the z-direction 
funFlux   = @F_euler3d;                           % Function handle for flux vector computation in each direction 
funEig    = @Eig_euler3d;                         % Function handle for computing eigenvalues of flux jacobian    
funInit   = @Q0;                                  % Function handle for initial conditions                        
funBC     = @BCperiodic;                          % Function handle for boundary conditions                       
isSource  = true;                                 % Flag to use user-defined "sourceTerm.m" to handle source terms    
gam       = 1.4;                                  % Ratio of specific heats
isConsErr = false;                                % Flag to calculate the local conservation error

%% For spead measurements, make sure to use only a single thread
maxNumCompThreads(1); 
fprintf("Using %d threads\n",maxNumCompThreads);

%% Setup library paths
%restoredefaultpath;
%run(getenv("TTWENOLIB")+"/ft_setup.m");
%% Setup equation specific path for flux and eigenvalue routines 
%addpath(fdwenolib.ft+"/equations/Euler3d");

%% Setup arrays for error and speed measurement reporting
L1error   = zeros(Neq,num_lvls);
L2error   = zeros(Neq,num_lvls);
Linferror = zeros(Neq,num_lvls);
CPUtime   = zeros(num_lvls,1);
N         = zeros(num_lvls,3); 

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
    Ny = Nx;
    Nz = Nx;
    %
    N(lvl,:) = [Nx Ny Nz];
    %
    fprintf("%s",repelem("-",70));fprintf("\n");
    fprintf("Starting to solve grid level=%d with Nx=%d Ny=%d Nz=%d\n",lvl,Nx,Ny,Nz);
    fprintf("%s",repelem("-",70));fprintf("\n");
    %
    data = ftFDWenoSolver(Tfinal, dtOption, CFL, Neq, Nx, Ny, Nz, Lx, Ly, Lz,...
                          funFlux, funEig, funInit, funBC, isSource, gam, isConsErr);
    % Set time spent by the solver
    CPUtime(lvl) = data.CPUtime;
    %% Error calculation
    Q  = data.Q;                                         % numerical result
    Qe = solExact(data.X,data.Y,data.Z,Tfinal,data.gam); % exact solution
    %
    vol = Lx*Ly*Lz/(Nx*Ny*Nz);
    h32 = sqrt(vol);
    %
    for eq=1:Neq
        %
        dQ = abs(Q(eq,:)-Qe(eq,:));
        %
        L1error(eq,lvl)   = sum(dQ)*vol;
        %
        L2error(eq,lvl)   = norm(dQ)*h32;
        %
        Linferror(eq,lvl) = max(dQ);
        %
    end
    %
    save_errors(resultPath+"/L1.txt",L1error(:,1:lvl));
    save_errors(resultPath+"/L2.txt",L2error(:,1:lvl));
    save_errors(resultPath+"/Linf.txt",Linferror(:,1:lvl));
    save_errors(resultPath+"/time.txt",CPUtime(1:lvl));
    %
    print_errors(1,L1error(:,1:lvl),"L1",N(1:lvl,:));
    print_errors(1,L2error(:,1:lvl),"L2",N(1:lvl,:));
    print_errors(1,Linferror(:,1:lvl),"Linf",N(1:lvl,:));
    %
    fprintf("\n");
end
%
save_error_table(resultPath+"/Table-L1.txt",L1error,"L1",N);
save_error_table(resultPath+"/Table-L2.txt", L2error,"L2",N);
save_error_table(resultPath+"/Table-Linf.txt",Linferror,"Linf",N);
%
fprintf("\n ------ ------------\n");
fprintf(" level | CPUtime (s)\n");
fprintf(" ------ ------------\n");
fprintf(" %5d | %10.4e\n",[(1:length(CPUtime))' CPUtime]');
fprintf("\n");