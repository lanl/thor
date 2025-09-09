function LUtt = tt_laplacian_3d(Utt, dx, dy, dz, diffcoeff, tt_tol)
%LaplacianUtt = tt_laplacian_3d(Utt, dx, dy, dz, diffcoeff, tt_tol)
% tensor version of laplacian_3d_dirichlet
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
Nx = Utt.n(1);
Ny = Utt.n(2);
Nz = Utt.n(3);

G = core2cell(Utt);
%FD on the 1st core
G{1}(:,2:Nx-1,:) = (G{1}(:,1:Nx-2,:) - 2*G{1}(:,2:Nx-1,:) + G{1}(:,3:Nx,:))/dx/dx;
Uxxtt = cell2core(tt_tensor,G);
% Uxxtt = tt_set_zero_boundaries(Uxxtt,{[1,Nx],[1,Ny],[1,Nz]});

% Uyytt
G = core2cell(Utt);
%FD on the 2nd core
G{2}(:,2:Ny-1,:) = (G{2}(:,1:Ny-2,:) - 2*G{2}(:,2:Ny-1,:) + G{2}(:,3:Ny,:))/dy/dy;
Uyytt = cell2core(tt_tensor,G);
% Uyytt = tt_set_zero_boundaries(Uyytt,{[1,Nx],[1,Ny],[1,Nz]});


% Uzztt
G = core2cell(Utt);
%FD on the 3rd core
G{3}(:,2:Nz-1,:) = (G{3}(:,1:Nz-2,:) - 2*G{3}(:,2:Nz-1,:) + G{3}(:,3:Nz,:))/dz/dz;
Uzztt = cell2core(tt_tensor,G);
% Uzztt = tt_set_zero_boundaries(Uzztt,{[1,Nx],[1,Ny],[1,Nz]});


% Laplacian U
LUtt = round(diffcoeff*(Uxxtt + Uyytt + Uzztt), tt_tol);

% set LU's boundary to zeroes (routine in tt_func)
LUtt = tt_set_zero_boundaries(LUtt,{[1,Nx],[1,Ny],[1,Nz]});

end