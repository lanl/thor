function dU = fg_oderhs(t,U,prob) 
% zero boundary conditions 
% Unpack problem info 
dx = prob.deltas(1); dy = prob.deltas(2); dz = prob.deltas(3); 
diffcoeff = prob.diffcoeff; 
fn_forcing = prob.forcingfn_xyz;
mask = prob.mask;

% Set boundary values for U 
U_ext = zeros(size(U));
U(mask) = U_ext(mask);

% Compute laplacian 
LaplacianU = laplacian_3d_dirichlet(U, dx, dy, dz, diffcoeff);

% Forcing term
Fterm = fn_forcing(t);

dU = LaplacianU + Fterm;

% Overwrite boundary
%    Here we set dU(i,j,k) = 0 on all boundary faces
dU_ext = zeros(size(U));
dU(mask) = dU_ext(mask);
end 