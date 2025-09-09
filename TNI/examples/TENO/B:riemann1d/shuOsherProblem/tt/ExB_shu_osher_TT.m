%% Info
% Example B   : Shu-Osher problem in tensor-train format
% Developer   : M. Engin Danis
% Notes       :
% It is recommended to run this script on terminal:
% (1) do: export TTTENOLIB=/path/of/TENO-paper
% (2) do: matlab -nodesktop -nojvm -nosplash -nodisplay -singleCompThread -r "try, run('ExB_shu_osher_TT.m'); catch err, disp(getReport(err)); end;exit;"
% This script writes 1d results (x, rho_tt, u_tt, p_tt) for each grid level under "results" directory, which can be used for postprocessing if needed.
addpath(genpath('../../../../../matlab/TENO/src/tt/'))
addpath(genpath('../../../../../matlab/TENO/src/misc/'))
addpath(genpath('../../../../../matlab/TENO/utils/'))
addpath(genpath('../../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../../matlab/utils/ttfunc/'))

%% set number of grid levels 
num_lvls = 4;
 
%% Input parameters to ttFDTenoSolver.m (Nx, Ny, Nz will be set later)
Tfinal    = 1.8;                                  % Final Time of the simulation   
dtOption  = 3;                                    % Time Stepping Option           
CFL       = 0.05;                                 % CFL Number for adaptive time-stepping
Neq       = 5;                                    % Number of equations 
N1d       = 50*(2.^linspace(1,num_lvls,num_lvls));% Number of cells in each direction
Lx        = 10;                                   % Domain length in the x-direction 
Ly        = 10;                                   % Domain length in the y-direction 
Lz        = 10;                                   % Domain length in the z-direction 
funFlux   = @F_euler3d;                           % Function handle for flux vector computation in each direction 
funEig    = @Eig_euler3d;                         % Function handle for computing eigenvalues of flux jacobian    
funInit   = @Q0;                                  % Function handle for initial conditions                        
funBC     = @BC;                                  % Function handle for boundary conditions                       
ReconChar = true;                                 % Flag to use characteristic decomposition (true or false) 
isSource  = false;                                % Flag to use user-defined "sourceTerm.m" to handle source terms    
gam       = 1.4;                                  % Ratio of specific heats (not used here, so it is left empty)
eps_tt    = 1e-14;                                % Initial tolerance for tt trunctation                                   
eps_cross = 1e-14;                                % Initial tolerance for tt cross interpolation                           
C_eps     = 2/sqrt(Lx*Ly*Lz);                     % Empricial coefficient in Eq. (35) of the TT-TENO paper by Danis et al. 
isConsErr = false;                                % Flag to calculate the local conservation error

%% For spead measurements, make sure to use only a single thread
maxNumCompThreads(1); 
fprintf("Using %d threads\n",maxNumCompThreads);

%% Setup arrays for error and speed measurement reporting
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
    Nz = Nx; 
    %
    fprintf("%s",repelem("-",70));fprintf("\n");
    fprintf("Starting to solve grid level=%d with Nx=%d Ny=%d Nz=%d\n",lvl,Nx,Ny,Nz);
    fprintf("%s",repelem("-",70));fprintf("\n");
    %
    data = ttFDTenoSolver(Tfinal, dtOption, CFL, Neq, Nx, Ny, Nz, Lx, Ly, Lz,...
                          funFlux, funEig, funInit, funBC, ReconChar, isSource, gam, ...
                          eps_tt, eps_cross, C_eps, isConsErr);
    % Set time spent by the solver
    CPUtime(lvl) = data.CPUtime;
    % Write data to file
    Q  = data.ft.Q; % numerical result (converted full-tensor form)
    %
    rho = Q(1,:,1,1);
    u   = Q(2,:,1,1)./rho;
    p   = (gam-1)*(Q(5,:,1,1) - 0.5*rho.*u.^2);
    %
    fileID = fopen(sprintf("%s/u-lvl=%d.txt",resultPath,lvl),'w');
    %
    for i=1:Nx
        %
        fprintf(fileID,'%.15e\t%.15e\t%.15e\t%.15e\n',data.X(i,1,1),rho(i),u(i),p(i));
        %
    end
    fclose(fileID);
    %
    save_errors(resultPath+"/time.txt",CPUtime(1:lvl));
    %
    save_ranks(resultPath+"/ranks-lvl="+lvl+".txt",data.tt.rank);
    %
    fprintf("\n");
end
%
fprintf("\n ------ ------------\n");
fprintf(" level | CPUtime (s)\n");
fprintf(" ------ ------------\n");
fprintf(" %5d | %10.4e\n",[(1:length(CPUtime))' CPUtime]');
fprintf("\n");