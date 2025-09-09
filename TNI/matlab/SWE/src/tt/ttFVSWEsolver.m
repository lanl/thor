% ---------------------------------------------------------------------------------------- %
% Developer   : Mustafa Engin Danis (@LANL, T-1)                                           %
% Date        : 05/01/2024                                                                 %
% Last Update : 07/01/2025                                                                 %
% Purpose     : Full-tensor Finite Volume Solver for Shallow Water Equations               %
% Publication : Danis, M. Engin et al.,                                                    %
%               High-order Tensor-Train Finite Volume Methods for Shallow Water Equations  %
%               AMS Monthly Weather Review                                                 %
%               https://doi.org/10.1175/MWR-D-24-0165.1                                    %
% Please cite above paper if you use any part of this library                              %
% ---------------------------------------------------------------------------------------- %
% Input     | Description                                                                  %
% ---------------------------------------------------------------------------------------- %
% SWEtype   : Shallow Water Equations (SWEs) type                                          %
%              SWEtype = "linear" uses linear SWEs                                         %
%              SWEtype = "nonlinear" uses nonlinear SWEs                                   %
% Tfinal    : Final Time of the simulation                                                 %
% dtOption  : Time Stepping Option                                                         %
%              1 = standard dt = cfl*dx_min/max_eig                                        %
%              2 = conservative dt = cfl*(dx_min^(5/3)) to maitain                         %
%                  5th order accuracy also in time                                         %
%              3 = simple dt = cfl*dx_min for debugging purposes                           %
% CFL       : CFL Number for adaptive time-stepping                                        %
% C_eps     : Empricial coefficient in Eq. (35) of the TT-WENO paper by Danis et al.  % 
% Nx        : Number of cells in the x-direction (interior)                                %
% Ny        : Number of cells in the y-direction (interior)                                %
% Lx        : Domain length in the x-direction                                             %
% Ly        : Domain length in the y-direction                                             %
% H         : Mean depth only used for linear SWEs                                         %
% gacc      : Acceleration of gravity                                                      %
% Coriolis_f: Coriolis parameter                                                           %
% Q0        : Function handle for initial conditions                                       %
% BC        : Function handle for boundary conditions                                      %
% Recon     : Reconstruction method type ("upwind3", "upwind5" or "weno5")                 %
% Ref       : Struct for reference values. Must include:                                   %
%              - Ref.L : reference length scale                                            %
%              - Ref.U : reference velocity scale                                          %
%              - Ref.H : reference height scale                                            %
% isManuf   : Flag to use user-defined "manuf.m" to handle manufactured solution           %
%             true  = use "manuf.m"                                                        %
%             false = do not use "manuf.m"                                                 %
%             "manuf.m": must be located in the directory of an example                    %
% isConsErr : Flag to calculate the conservation error (true or false)                     %
% ---------------------------------------------------------------------------------------- %
% Output    |  Description                                                                 %
% ---------------------------------------------------------------------------------------- %
% data      : struct that contains several solver information including the results        %
% ---------------------------------------------------------------------------------------- %
function data = ttFVSWEsolver(SWEtype, Tfinal, dtOption, CFL, C_eps, Nx, Ny, Lx, Ly, H, ...
                              gacc, Coriolis_f, Q0, BC, Recon, Ref, isManuf, isConsErr)
    %
    % Function Structure:
    % 0. Creating the main data structure
    % 1. Initializing the meta data and mesh
    % 2. Allocating memory  
    % 3. Initiazlizing the solution vector
    % 4. Timestepping
    % 5. Preparing the solution to return the output 
    %
    % Creating the main data structure
    data = {};
    % allow passing an empty argument for the nonlinear case
    if(SWEtype == "nonlinear")
        if(~exist("H","var"))
            H = 0; 
        elseif(isempty(H))
            H = 0;
        end
    end
    %
    % Initializing the solver data
    %
    data = initMetaData(data, SWEtype, Tfinal, dtOption, CFL, C_eps, Nx, Ny, Lx, Ly, H,...
                        gacc, Coriolis_f, Q0, BC, Recon, Ref, isManuf, isConsErr);
    %
    % Allocating memory
    %
    data = allocateMem(data);
    %
    % Initialize the solution
    %
    data = data.Q0(data);
    %
    data = calculateEps(data);
    %
    data.tt.Q = applyBC(data.tt.Q,data,data.curr_time);
    %
    % Advance solution to time Tfinal
    %
    tic;
    data = explicitSolve(data);
    data.CPUtime=toc;
    %
    % Prepare the output
    %
    data = prepareOutput(data);
    %
end