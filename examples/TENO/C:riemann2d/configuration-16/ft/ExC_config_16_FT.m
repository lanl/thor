%% Info
% Example C   : 2D Riemann problem (Configuration 16) in full-tensor format
% Developer   : M. Engin Danis
% Notes       :
% It is recommended to run this script on terminal:
% (1) do: export TTTENOLIB=/path/of/TENO-paper
% (2) do: matlab -nodesktop -nojvm -nosplash -nodisplay -singleCompThread -r "try, run('ExC_config_16_FT.m'); catch err, disp(getReport(err)); end;exit;"
% This script writes 2d results (x, y, rho_ft, rhou_ft, rhov_ft, rhoE_ft) for each grid level under "results" directory, which can be used for postprocessing if needed.
addpath(genpath('../../../../../matlab/TENO/src/ft/'))
addpath(genpath('../../../../../matlab/TENO/src/misc/'))
addpath(genpath('../../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../../matlab/utils/ttfunc/'))

%% set number of grid levels 
num_lvls = 1;
 
%% Input parameters to ftFDTenoSolver.m (Nx, Ny, Nz will be set later)
Tfinal    = 0.2;                                   % Final Time of the simulation   
dtOption  = 1;                                     % Time Stepping Option           
CFL       = 0.475;                                 % CFL Number for adaptive time-stepping
Neq       = 5;                                     % Number of equations 
N1d       = 200*(2.^linspace(1,num_lvls,num_lvls));% Number of cells in each direction
Lx        = 1;                                     % Domain length in the x-direction 
Ly        = 1;                                     % Domain length in the y-direction 
Lz        = 1;                                     % Domain length in the z-direction 
funFlux   = @F_euler3d;                            % Function handle for flux vector computation in each direction 
funEig    = @Eig_euler3d;                          % Function handle for computing eigenvalues of flux jacobian    
funInit   = @Q0;                                   % Function handle for initial conditions                        
funBC     = @BC;                                   % Function handle for boundary conditions                       
ReconChar = true;                                  % Flag to use characteristic decomposition (true or false) 
isSource  = false;                                 % Flag to use user-defined "sourceTerm.m" to handle source terms    
gam       = 1.4;                                   % Ratio of specific heats (not used here, so it is left empty)
isConsErr = false;                                 % Flag to calculate the local conservation error

%% For spead measurements, make sure to use only a single thread
maxNumCompThreads(1); 
fprintf("Using %d threads\n",maxNumCompThreads);


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
    Ny = Nx; 
    Nz = 1;  % change this to Nx for a full 3d simulation
    %
    fprintf("%s",repelem("-",70));fprintf("\n");
    fprintf("Starting to solve grid level=%d with Nx=%d Ny=%d Nz=%d\n",lvl,Nx,Ny,Nz);
    fprintf("%s",repelem("-",70));fprintf("\n");
    %
    data = ftFDTenoSolver(Tfinal, dtOption, CFL, Neq, Nx, Ny, Nz, Lx, Ly, Lz,...
                          funFlux, funEig, funInit, funBC, ReconChar, isSource, gam, isConsErr);
    % Set time spent by the solver
    CPUtime(lvl) = data.CPUtime;
    % Write data to file
    fileID = fopen(sprintf("%s/Q-lvl=%d.txt",resultPath,lvl),'w');
    %
    for j=1:Ny
        for i=1:Nx
            %
            fprintf(fileID,'%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n',data.X(i,j,1),data.Y(i,j,1),data.Q(1,i,j,1),data.Q(2,i,j,1),data.Q(3,i,j,1),data.Q(5,i,j,1));
            %
        end
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