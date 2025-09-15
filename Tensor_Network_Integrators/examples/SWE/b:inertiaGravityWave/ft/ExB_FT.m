%% Info
% Example B   : Inertia-Gravity Wave in full-tensor format
% Developer   : M. Engin Danis
% Notes       :
% It is recommended to run this script on terminal:
% (1) do: export TTSWELIB=/path/of/SWE-paper
% (2) do: matlab -nodesktop -nojvm -nosplash -nodisplay -singleCompThread -r "try, run('ExB_FT.m'); catch err, disp(getReport(err)); end;exit;"

addpath(genpath('../../../../matlab/SWE/src/ft/'))
addpath(genpath('../../../../matlab/SWE/src/misc/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))

clear; clc; close all;

%% set the list of reconstruction schemes
ReconList = {"Upwind3", "Upwind5", "WENO5"};

%% set number of grid levels 
num_lvls = 5;

%% Input parameters to ftFVSWEsolver.m (Nx, Ny, CFL, Recon will be set later)
SWEtype    = "linear";                              % Shallow Water Equations (SWEs) type
Tfinal     = 10800;                                 % Final Time of the simulation
dtOption   = 2;                                     % Time Stepping Option
N1d        = 10*(2.^linspace(1,num_lvls,num_lvls)); % Number of cells in each direction
Lx         = 1e7;                                   % Domain length in the x-direction
Ly         = 1e7;                                   % Domain length in the y-direction
H          = 1000;                                  % Mean depth
gacc       = 10;                                    % Acceleration of gravity
Coriolis_f = 1e-4;                                  % Coriolis parameter
funInit    = @Q0;                                   % Function handle for initial conditions
funBC      = @BC;                                   % Function handle for boundary conditions 
Ref.L      = 1e7;                                   % Reference length scale
Ref.H      = 0.2;                                   % Reference depth scale
isManuf    = false;                                 % Flag to use user-defined "manuf.m" to handle manufactured solutions terms
isConsErr  = false;                                 % Flag to calculate conservation errors

%% Reference velocity scale is calculated as below
eta   = 0.2;
kx    = 4e-7*pi;
C2    = 1e4;
k     = sqrt(2)*kx;
om    = sqrt(C2*k^2+Coriolis_f^2);
Ref.U = gacc*eta*kx*om/(om^2-Coriolis_f^2);

%% For spead measurements, make sure to use only a single thread
maxNumCompThreads(1); 
fprintf("Using %d threads\n",maxNumCompThreads);

% Set number of equations
Neq = 3;

%% Setup arrays for error and speed measurement reporting
L1error   = zeros(Neq,num_lvls);
L2error   = zeros(Neq,num_lvls);
Linferror = zeros(Neq,num_lvls);
CPUtime   = zeros(num_lvls,1);
N         = zeros(num_lvls,2);

%% Set the path where results will be written. This directory will be created if needed
resultPath = "results";
if exist(resultPath, 'dir')
    cmd = "rm -rf "+resultPath;
    system(cmd);
end
mkdir(resultPath);

%% Run simulations for all reconstruction methods
for recon_idx = 1:length(ReconList)
    %% Get the reconstruction method type
    Recon = ReconList{recon_idx};

    %% Set the path where specific results for this reconstruction will be written 
    resultPathRecon = resultPath + "/" + Recon;
    if exist(resultPathRecon, 'dir')
        cmd = "rm -rf "+resultPathRecon;
        system(cmd);
    end
    mkdir(resultPathRecon);

    %% Set the CFL number for this reconstruction method
    if(Recon=="Upwind3")
        CFL = 0.0001;
    else
        CFL = 0.001;
    end

    %% Run simulations for all grid levels
    for lvl=1:length(N1d)
        %
        Nx = N1d(lvl);
        Ny = Nx;
        %
        N(lvl,:) = [Nx Ny];
        %
        fprintf("%s",repelem("-",70));fprintf("\n");
        fprintf("Starting to solve grid level=%d with Nx=%d Ny=%d\n",lvl,Nx,Ny);
        fprintf("%s",repelem("-",70));fprintf("\n");
        %
        data = ftFVSWEsolver(SWEtype, Tfinal, dtOption, CFL, Nx, Ny, Lx, Ly, H, gacc, ...
                             Coriolis_f, funInit, funBC, Recon, Ref, isManuf, isConsErr);
        % Set time spent by the solver
        CPUtime(lvl) = data.CPUtime;
        %% Error calculation
        Q  = data.Q;                                                      % numerical result
        Qe = applyQuadrature(@solExact,data,Tfinal/data.ref.t)/data.vol;  % exact solution
        %
        % Dimensionalize Qe
        %
        Qe(1,:,:) = Qe(1,:,:) * data.Qref(1);
        Qe(2,:,:) = Qe(2,:,:) * data.Qref(2);
        Qe(3,:,:) = Qe(3,:,:) * data.Qref(3);
        %
        vol = Lx*Ly/(Nx*Ny);
        h32 = sqrt(vol);
        %
        % Note: L1 and L2 error values might seem large due to very large Lx and Ly values.
        %       However, the convergence rate will recover the formal order of the underlying scheme.
        %       So, scaling them by Lx*Ly will produce more familiar error value (optional). 
        %       However, for testing purposes, L1 and L2 errors are calculated as below.
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
        save_errors(resultPathRecon+"/L1.txt",L1error(:,1:lvl));
        save_errors(resultPathRecon+"/L2.txt",L2error(:,1:lvl));
        save_errors(resultPathRecon+"/Linf.txt",Linferror(:,1:lvl));
        save_errors(resultPathRecon+"/time.txt",CPUtime(1:lvl));
        %
        print_errors(1,L1error(:,1:lvl),"L1",N(1:lvl,:));
        print_errors(1,L2error(:,1:lvl),"L2",N(1:lvl,:));
        print_errors(1,Linferror(:,1:lvl),"Linf",N(1:lvl,:));
        %
        fprintf("\n");
    end
    %
    save_error_table(resultPathRecon+"/Table-L1.txt",L1error,"L1",N);
    save_error_table(resultPathRecon+"/Table-L2.txt", L2error,"L2",N);
    save_error_table(resultPathRecon+"/Table-Linf.txt",Linferror,"Linf",N);
    %
    fprintf("\n ------ ------------\n");
    fprintf(" level | CPUtime (s)\n");
    fprintf(" ------ ------------\n");
    fprintf(" %5d | %10.4e\n",[(1:length(CPUtime))' CPUtime]');
    fprintf("\n");
    %
end

