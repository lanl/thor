% driver for heat spatial scaling problem 
% Run FG and TT

% problem inputs
% tspan - time interval
% nvals - vector of spatial grid sizes
% rtol - relative accuracy tolerance 
% atol - absolute tolerance
% tt_tol - rounding tolerance 
% mname - name of explicit Runge Kutta method
% uexact - anonymous function 

nvals = 40:40:120;

input.tspan = [0,0.2]; 
input.nvals = nvals;
input.diff = 0.01;                         % diffusion coefficient 
input.rtol = 1e-5; 
input.atol = 1e-12;
input.tt_tol = input.rtol; 
input.mname = 'Bogacki-Shampine-ERK';

input.uexact = @(t,x,y,z) sin(pi*x).* sin(pi*x).* sin(pi*y).* sin(pi*y).* ...
    sin(pi*z).* sin(pi*z).* cos(pi*t).* cos(pi*t);

heat_spatial_scaling_struct(input,['h3D_driver_',input.mname],1)