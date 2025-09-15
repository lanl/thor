function [y,lits,ierr] = newton(Fcn, Afn, y0, Fdata, ftol, stol, maxit)
% usage: [y,lits,ierr] = newton(Fcn, Afn, y0, Fdata, ftol, stol, maxit)
%
% Newton solver for the root-finding problem defined by the function Fcn,
%     F(y,Fdata) = 0
%
% Inputs:  Fcn = function name for nonlinear residual, F(y,Fdata).  Note, all
%                data required to evaluate F (other than y) should be stored
%                in the data structure Fdata.
%          Afn = function name for Jacobian of nonlinear residual,
%                A = partial_y F(y,Fdata).
%                Afn should use the same data structure for additional data
%                as F.
%          y0 = initial guess
%          Fdata = structure containing extra information for evaluating F.
%          ftol = desired nonlinear residual tolerance, ||F(y)|| < ftol
%          stol = desired solution tolerance, ||s|| < stol
%          maxit = maximum allowed iterations
% Outputs: y = solution to root-finding problem
%          lits = total # of linear solves taken
%          ierr = output flag denoting success (0) or failure (1)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% check solver inputs
if (maxit < 1)
   error('newton error: requires at least 1 iteration (maxit)');
end
if (stol <= 0)
   error('newton error: tolerance must be positive (stol)');
end
if (ftol <= 0)
   error('newton error: tolerance must be positive (ftol)');
end

% initialize result, increment vector, residual, statistics
y = y0;
s = ones(size(y));
lits = 0;

% perform iterations
for i=1:maxit

   % compute residual at current guess
   F = Fcn(y,Fdata);

   % check residual and increment for stopping
   if ((norm(s,inf) < stol) | (norm(F,inf) < ftol))
      ierr = 0;
      return
   end

   % compute Jacobian
   A = Afn(y,Fdata);

   % perform Newton update
   % [s,~] = gmres(A,F,10,ftol);
   % [s,~] = pcg(A,F,ftol);
   % [s,~] = bicgstab(A,F,ftol);
   s = A\F;
   y = y - s;
   lits = lits + 1;

end

% if we've made it to this point, the Newton iteration did not converge
ierr = 1;
%fprintf('\nnewton warning: nonconvergence after %i iterations (|F| = %g)\n',maxit,norm(F,inf));

% end of function
