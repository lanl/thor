function dUtt = tt_oderhs(t,Utt,prob,tt_tol) 
% Unpack problem info 
dx = prob.deltas(1); dy = prob.deltas(2); dz = prob.deltas(3); 
diffcoeff = prob.diffcoeff; 
Ctt = prob.Ctt;
boundaryfn = @(x,y,z) prob.boundaryfn(t,x,y,z);
fn_forcingtt = @(x,y,z) prob.forcingfn(t,x,y,z);

Nx = Utt.n(1); 
Ny = Utt.n(2); 
Nz = Utt.n(3); 

% Set boundary values for Utt 
Utt = tt_set_zero_boundaries(Utt,{[1,Nx],[1,Ny],[1,Nz]});

% Compute laplacian 
LaplacianUtt = tt_laplacian_3d(Utt, dx, dy, dz, diffcoeff,tt_tol);

% Forcing term 
Ftermtt = amen_cross_zero(Ctt, @(x) cross_fun_nD(x,fn_forcingtt),tt_tol,'verb',0);
Ftermtt = tt_set_zero_boundaries(Ftermtt,{[1,Nx],[1,Ny],[1,Nz]});

% Combine right hand side function 
dUtt = LaplacianUtt + Ftermtt;
dUtt = tt_set_zero_boundaries(dUtt,{[1,Nx],[1,Ny],[1,Nz]});

end 