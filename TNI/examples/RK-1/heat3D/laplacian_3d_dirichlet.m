function LaplacianU = laplacian_3d_dirichlet(U, dx, dy, dz, diffcoeff)
%LaplacianU = laplacian_3d_dirichlet(U, dx, dy, dz, diffcoeff)
% Computes the 3D Laplacian with Dirichlet boundary conditions
% Boundary values are part of U 
%
% INPUTS:
% U: Nx x Ny x Nz matrix (grid values, including boundaries)
% dx: Grid spacing in x-direction
% dy: Grid spacing in y-direction
% dz: Grid spacing in z-direction
% diffcoeff: Diffusion coefficient
%
% OUTPUTS:
% 3D LaplacianU : Nx x Ny x Nz matrix (same size as U, boundaries set to exact solution)

% Get grid dimensions
[Nx, Ny, Nz] = size(U);

% Initialize Laplacian
LaplacianU = zeros(Nx, Ny, Nz);

% Interior indices for finite difference computation
k = 2:Nz-1; % z-direction
j = 2:Ny-1; % y-direction
i = 2:Nx-1; % x-direction

% Compute Laplacian for interior points
LaplacianU(i, j, k) = diffcoeff * (...
    (U(i+1, j, k) - 2*U(i, j, k) + U(i-1, j, k)) / dx^2 + ... % x-direction
    (U(i, j+1, k) - 2*U(i, j, k) + U(i, j-1, k)) / dy^2 + ... % y-direction
    (U(i, j, k+1) - 2*U(i, j, k) + U(i, j, k-1)) / dz^2);     % z-direction

% LaplacianU is already initialized to zero so boundary values remain zero.
end