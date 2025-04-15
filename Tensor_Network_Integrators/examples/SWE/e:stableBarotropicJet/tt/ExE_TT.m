%% Info
% Example E   : Stable Barotropic Jet in tensor-train format
% Developer   : M. Engin Danis
% Notes       :
% It is recommended to run this script on terminal:
% (1) do: export TTSWELIB=/path/of/SWE-paper
% (2) do: matlab -nodesktop -nojvm -nosplash -nodisplay -singleCompThread -r "try, run('ExE_TT.m'); catch err, disp(getReport(err)); end;exit;"

addpath(genpath('../../../../matlab/SWE/src/tt/'))
addpath(genpath('../../../../matlab/SWE/src/misc/'))
addpath(genpath('../../../../matlab/SWE/utils/'))
addpath(genpath('../../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../../matlab/utils/ttfunc/'))

clear; clc; close all;

%% set the list of reconstruction schemes
ReconList = {"Upwind3", "Upwind5", "WENO5"};

%% set number of grid levels 
num_lvls = 7;

%% Input parameters to ttFVSWEsolver.m (Nx, Ny, CFL, C_eps, Recon will be set later)
SWEtype    = "nonlinear";                           % Shallow Water Equations (SWEs) type
Tfinal     = 3600*24*5;                             % Final Time of the simulation
dtOption   = 2;                                     % Time Stepping Option
N1d        = 20*(2.^linspace(1,num_lvls,num_lvls)); % Number of cells in each direction
Lx         = 2*pi*6371220;                          % Domain length in the x-direction
Ly         = Lx/2;                                  % Domain length in the y-direction
H          = [];                                    % Mean depth (not needed in this example)
gacc       = 9.80616;                               % Acceleration of gravity
Coriolis_f = 2*(7.292e-5)*sin(pi/4);                % Coriolis parameter
funInit    = @Q0;                                   % Function handle for initial conditions
funBC      = @BC;                                   % Function handle for boundary conditions 
Ref.L      = Lx;                                    % Reference length scale
Ref.U      = 20;                                    % Reference velocity scale
Ref.H      = 8000;                                  % Reference depth scale
isManuf    = false;                                 % Flag to use user-defined "manuf.m" to handle manufactured solutions terms
isConsErr  = true;                                  % Flag to calculate conservation errors

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
        CFL   = 0.005;
        C_eps = 1;
    elseif(Recon=="Upwind5")
        CFL   = 0.5;
        C_eps = 1;
    else
        CFL   = 0.5;
        C_eps = 1;
    end

    %% Run simulations for all grid levels
    for lvl=1:length(N1d)
        %
        Nx = N1d(lvl);
        Ny = Nx/2;
        %
        N(lvl,:) = [Nx Ny];
        %
        fprintf("%s",repelem("-",70));fprintf("\n");
        fprintf("Starting to solve grid level=%d with Nx=%d Ny=%d\n",lvl,Nx,Ny);
        fprintf("%s",repelem("-",70));fprintf("\n");
        %
        data = ttFVSWEsolver(SWEtype, Tfinal, dtOption, CFL, C_eps, Nx, Ny, Lx, Ly, H, gacc, ...
                             Coriolis_f, funInit, funBC, Recon, Ref, isManuf, isConsErr);
        % Set time spent by the solver
        CPUtime(lvl) = data.CPUtime;
        %% Error calculation
        Q  = data.Q;       % numerical result
        Qe = data.tt.Q0;   % exact solution (this is actually the same as the initial condition)
        %
        Qe = remove_ghost(Qe,Neq);
        %
        % Dimensionalize Qe
        %
        Qe{1} = Qe{1} * data.Qref(1);
        Qe{2} = Qe{2} * data.Qref(2);
        Qe{3} = Qe{3} * data.Qref(3);
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
            dQ = reshape(abs(full(Q{eq})-full(Qe{eq})),[Nx*Ny,1]);
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
        save_ranks(resultPathRecon+"/ranks-lvl="+lvl+".txt",data.tt.rank);
        %
        print_errors(1,L1error(:,1:lvl),"L1",N(1:lvl,:));
        print_errors(1,L2error(:,1:lvl),"L2",N(1:lvl,:));
        print_errors(1,Linferror(:,1:lvl),"Linf",N(1:lvl,:));
        %
        if(isConsErr)
            %
            save_errors(sprintf("%s/conservation-int-lvl=%d.txt",resultPathRecon,lvl),data.consInt);
            % normalize as well
            data.consInt(:,2) = (data.consInt(:,2)-data.consInt(1,2))/data.consInt(1,2);
            data.consInt(:,3) = (data.consInt(:,3)-data.consInt(1,3))/data.consInt(1,3);
            %
            save_errors(sprintf("%s/conservation-err-lvl=%d.txt",resultPathRecon,lvl),data.consInt);
            %
        end
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

