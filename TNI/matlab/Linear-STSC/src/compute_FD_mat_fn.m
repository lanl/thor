function [A,Tvec] = compute_FD_mat_fn(N,a,b)

dx = (b-a)/(N-1);
A = 1/dx*(diag(1*ones(N,1)) + diag(-1*ones(N-1,1),-1));
Tvec = a:dx:b;