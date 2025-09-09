% ---------------------------------------------------------------------------------- %
% Developer   : Mustafa Engin Danis (@LANL, T-1)                                     %
% Date        : 10/26/2023                                                           %
% Last Update : 04/15/2025                                                           %
% Purpose     : Full-tensor Finite Difference Weighted Essentially                   %
%               Non-oscillatory (WENO) Solver for Hyperbolic PDEs                    %
% Publication : Danis, M. Engin et al.,                                              %
%               Tensor-Train WENO Scheme for Compressible Flows,                     %
%               Journal of Computational Physics 529.C (2025)                        %
%               https://doi.org/10.1016/j.jcp.2025.113891                            %
% Please cite above paper if you use any part of this library                        %
% ---------------------------------------------------------------------------------- %
% Input    | Description                                                             %
% ---------------------------------------------------------------------------------- %
% Tfinal   : Final Time of the simulation                                            %
% dtOption : Time Stepping Option                                                    %
%             1 = standard dt = cfl*dx_min/max_eig                                   %
%             2 = conservative dt = cfl*(dx_min^(5/3)) to maitain                    %
%                 5th order accuracy also in time                                    %
%             3 = simple dt = cfl*dx_min for debugging purposes                      %
% CFL      : CFL Number for adaptive time-stepping                                   %
% Neq      : Number of equations (e.g. Neq=5 for 3D Euler )                          % 
% Nx       : Number of cells in the x-direction (interior)                           %
% Ny       : Number of cells in the y-direction (interior)                           %
% Nz       : Number of cells in the z-direction (interior)                           %
% Lx       : Domain length in the x-direction                                        %
% Ly       : Domain length in the y-direction                                        %
% Lz       : Domain length in the z-direction                                        %
% F        : Function handle for flux vector computation in each direction           %
% Eig      : Function handle for computing eigenvalues of flux jacobian              %
% Q0       : Function handle for initial conditions                                  %
% BC       : Function handle for boundary conditions                                 %
% isSource : Flag to use user-defined "sourceTerm.m" to handle source terms          %
%            (e.g. manufactured solutions, buoyancy, etc.)                           %
%             true  = use "sourceTerm.m"                                             %
%             false = do not use "sourceTerm.m"                                      %
%             "sourceTerm.m": must be located in the directory of an example         %
% gam      : Ratio of specific heats                                                 %
% isConsErr: Flag to calculate the local conservation error (true or false)          %
% ---------------------------------------------------------------------------------- %
% Output |  Description                                                              %
% ---------------------------------------------------------------------------------- %
% data   : struct that contains several solver information including the results     %
% ---------------------------------------------------------------------------------- %
function data = ftFDWenoSolver(Tfinal, dtOption, CFL, Neq, Nx, Ny, Nz, Lx, Ly, Lz,...
                               F, Eig, Q0, BC, isSource, gam, isConsErr)
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
    %
    % Initializing the solver data
    %
    data = initMetaData(data, Tfinal, dtOption, CFL, Neq, Nx, Ny, Nz, Lx, Ly, Lz,...
                        F, Eig, Q0, BC, isSource, gam, isConsErr);
    %
    % Allocating memory
    %
    data = allocateMem(data);
    %
    % Initialize the solution
    %
    data = data.Q0(data);
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
