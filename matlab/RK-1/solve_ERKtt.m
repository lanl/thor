function [tvals,Ytt,nsteps,output] = solve_ERKtt(fcn,tvals,Y0,par)
% usage: [tvals,Ytt,nsteps,output] = solve_ERKtt(fcn,tvals,Y0,par)
%
% Adaptive time step explicit Runge-Kutta solver for the
% vector-valued ODE problem
%     y' = F(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     fcn    = string holding function name for F(t,Y)
%     tvals  = [t0, t1, t2, ..., tN]
%     Y0     = initial value tensor 
%     par    = struct with solver parameters 
%     par fields 
%     ----------------
%     B      = Butcher matrix for ERK coefficients, of the form
%                 B = [c A;
%                      q b;
%                      p b2 ]
%              Here, c is a vector of stage time fractions (s-by-1),
%                    A is a matrix of Butcher coefficients (s-by-s),
%                    q is an integer denoting the method order of accuracy,
%                    b is a vector of solution weights (1-by-s),
%                    p is an integer denoting the embedding order of accuracy,
%                    b2 is a vector of embedding weights (1-by-s),
%              The [p, b2] row is optional.  If that row is not
%              provided the method will default to taking fixed
%              step sizes of size dtmin.
%     rtol   = desired relative error of solution  (scalar)
%     atol   = desired absolute error of solution  (vector or scalar)
%     dtmin   = minimum internal time step size (dtmin <= t(i)-t(i-1), for all i)
%     dtmax   = maximum internal time step size (dtmax >= dtmin)
%     dtinit  = initial internal time step size (dtmin <= dtinit <= dtmax)
%     ettinit = initial truncation tolerance 
%
% Outputs:
%     tvals  = the same as the input array tvals
%     y      = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%               y(t*) is a column vector of length m.
%     nsteps = number of internal time steps taken by method
%     output = struct with solver stats
%     output fields 
%     -----------------
%     dtvals = time steps taken 
%     nfail = total number of step size failures 
%     ettvals = truncation tolerances used 
%
% Note: to run in fixed-step mode, call with dtmin=dtmax as the desired
% time step size.
%
% (please acknowledge in any additional iterations) 
%  Original code by Daniel Reynolds (https://github.com/drreynolds/rklab)
% 
% This version: Rujeko Chinomona

%% Unpack solver parameters
struct2vars(par);

%% Adaptivity setup
% determine whether adaptivity is desired
adaptive = 0;
if (abs(dtmax-dtmin)/abs(dtmax) > sqrt(eps))
   adaptive = 1;
end

% if adaptivity enabled, determine approach for error estimation,
% and set the lower-order of accuracy accordingly
[Brows, Bcols] = size(B);
embedded = 0;
p = 0;
if (dtmax > dtmin)      % check whether adaptivity is desired
   if (Brows > Bcols)
      if (max(abs(B(Bcols+1,2:Bcols))) > eps)   % check for embedding coeffs
         embedded = 1;
         p = B(Bcols+1,1);
      end
   end
end

%% Solver parameters
dt_reduce = 0.1;          % failed step reduction factor
dt_safety = 0.9;          % adaptivity safety factor
dt_growth = 10;           % adaptivity growth bound
ONEMSM   = 1-sqrt(eps);  % coefficients to account for
ONEPSM   = 1+sqrt(eps);  %   floating-point roundoff
ERRTOL   = 1.1;          % upper bound on allowed step error
                           %   (in WRMS norm)

%% Initialization
N = length(tvals);
Ytt = cell(N,1); % outputs
Ytt{1} = Y0;
% initialize diagnostics
a_fails = 0;   % total accuracy failures
% initialize temporary variables
t = tvals(1);
Ynew = Y0;
% set initial time step size
dt = dtinit;

% initialize work counter
nsteps = 0;

err_rel = [];
dtstore = [];
ettstore = [];
est = [];
err_abs = [];
ett = ettinit;

%% Iterate over output time steps
for tstep = 2:length(tvals)
   % loop over internal time steps to get to desired output time
   while ((t-tvals(tstep))*dt < 0)

      % bound internal time step
      dt = max([dt, dtmin]);            % enforce minimum time step size
      dt = min([dt, dtmax]);            % maximum time step size
      dt = min([dt, tvals(tstep)-t]);  % stop at output time
      dtstore = [dtstore,dt];
      ettstore = [ettstore,ett];

      % reset step failure flag
      st_fail = 0;

      % compute updated solution and error estimate (if possible);
      % increment internal time steps counter
      if (adaptive)
            [Ynew,Yerr] = tt_ERKstep_embedded(fcn, Y0, t, dt, B,ett);
            nsteps = nsteps + 1;
      else
         [Ynew] = tt_ERKstep_basic(fcn, Y0, t, dt, B,ett);
         nsteps = nsteps + 1;
      end

      % if time step adaptivity enabled, check step accuracy
      if (adaptive)

         % estimate error in current step
         err_step = max(norm(Yerr)/norm(rtol*Ynew + atol), eps);
         abs_err = norm(Yerr);
         rel_err = abs_err/norm(Ynew);

         % Store 
         err_rel = [err_rel,rel_err];
         est = [est, err_step];
         err_abs = [err_abs,abs_err];

         % if error too high, flag step as a failure (will be recomputed)
         if (err_step > ERRTOL*ONEPSM)
            a_fails = a_fails + 1;
            st_fail = 1;
         end

      end

      % if step was successful (i.e. error acceptable)
      if (st_fail == 0)

         % update solution and time for last successful step
         Y0 = Ynew;
         t  = t + dt;

         % for adaptive methods, use error estimate to adapt the time step
         if (adaptive)

            h_old = dt;
            if (err_step == 0.0)     % no error, set max possible
               dt = tvals(end)-t;
            else                     % set next dt (I-controller)
               dt = dt_safety * h_old * err_step^(-1.0/p);
            end

            % enforce maximum growth rate on step sizes
            dt = min(dt_growth*h_old, dt);

         % otherwise, just use the fixed minimum input step size
         else
            dt = dtmin;
         end
         
         % Experiment with adaptive ett (To come)

     % if error test failed
     else

        % if already at minimum step, just return with failure
        if (dt <= dtmin)
           error('Cannot achieve desired accuracy.\n  Consider reducing dtmin or increasing rtol.\n');
           return
        end

        % otherwise, reset guess, reduce time step, retry solve
        Ynew = Y0;
        dt    = dt * dt_reduce;

        % Experiment with adaptive ett (To come)
   
     end  % end logic tests for step success/failure
   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   
   Ytt{tstep} = Ynew;

end  % time step loop


output = struct('dtvals',dtstore,'nfail',a_fails,'ettvals',ettstore);

% figure;
% semilogy(1:length(dtstore),dtstore)
% hold on
% semilogy(1:length(err_rel),err_rel,'r-')
% hold on
% semilogy(1:length(est),est,'m-')
% hold on
% yline(rtol,'k')
% hold on
% yline(ERRTOL*ONEPSM,'k')
% hold on 
% semilogy(1:length(ettstore),ettstore,'g')

end % end solve_ERK function












