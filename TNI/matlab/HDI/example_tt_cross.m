% this cript will test some methods to approximate a scalar function on a
% multi-D grid using tensor train
% some methods available in the toolbox: TT-Toolbox-master/cross/
% (1) multifuncrs : CAR tt-cross??? ; (2) dmrg; (3) amen;
% greedy2_cross.m might implement another method.

% Questions to investigate: 
% (1) Which methods work the best with quadrature rules? Do we need to rewrite the
% function to take non-equispaced grid as an input?
% (2) Which method is the easiest to translate to multiprecision?

close all; clear; clc;
run setup.m;
% run setup_advanpix_drac.m;
run setup_advanpix_temp.m;

%% test dmrg methods

n=10; d=3; q=10;
% fun = @(ind) mp(sum(ind));
fun= @(x) mp((q+1)/2.*sum(x.^q));
xgrid = sort(mp(rand(1,n)));

tt= amen_cross_mp([n,n,n],fun,mp(1e-20),xgrid); %

%COMMENTS: % to take Gaussian nodes as input, this function needs to be
%rewrite
display(tt) %print out the form of the tensor train
f_tt = full(tt); %this returns a vector

% --- check the error
% (1) Evaluate the function on a full 3D grid 50x50x50 -- feval
feval = mp(zeros(n,n,n));
for i = 1:n
  for j = 1:n
    for k = 1:n
      feval(i,j,k) = fun(xgrid([i,j,k]));
%       feval(i,j,k) = fun([i,j,k]);
    end
  end
end

% (2) compare the error between tt and feval
Error = norm(full(f_tt)-feval(:)) % not good err;


%% test cross2d_new -- this can only be used for 2D cases.
n = 500;
m = 500;
f = @(i,j) 1.0./(i + 2 * j);
%f = @(i,j) 1.0;
[u, v] = cross2d_new(f, n, m, 1e-12);
f_tt = u*v';

%test the error
[X1,X2] = ndgrid([1:500]);
fprintf('Approximation error = %.5e \n',norm(f(X1,X2)-f_tt));

%COMMENTS: % to take Gaussian nodes as input, this function