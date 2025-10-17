function y = RK45(odefun, tspan, y0, h, tol, eps_tt)
  % Initialize variables
  t0 = tspan(1);
  tf = tspan(2);
  t  = t0;
  y  = y0;
  %
  % Butcher Tableau for RK45 (Dormand-Prince)
  %
  A = [0, 0, 0, 0, 0, 0;
       1/5, 0, 0, 0, 0, 0;
       3/40, 9/40, 0, 0, 0, 0;
       44/45, -56/15, 32/9, 0, 0, 0;
       19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0;
       9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0];
  %
  b4 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84];                    % 4th order
  b5 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];  % 5th order
  c  = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
  % 
  tstep = 1;
  %fprintf("t(%d)=%e\n",tstep,t);
  while (t < tf)
    %
    tstep = tstep + 1;
    %
    tn = t;
    %
    yn = y;
    % 
    if (tn + h > tf)
      h = tf - tn;
    end
    %
    k1 = odefun(tn, yn);
    %
    ys = round(yn + h * A(2,1) * k1, eps_tt);
    k2 = odefun(tn + c(2) * h, ys);
    %
    ys = round(yn + h * round(A(3,1) * k1 + A(3,2) * k2, eps_tt), eps_tt);
    k3 = odefun(tn + c(3) * h, ys);
    %
    ys = round(yn + h * round(A(4,1) * k1 + round(A(4,2) * k2 + A(4,3) * k3, eps_tt), eps_tt), eps_tt);
    k4 = odefun(tn + c(4) * h, ys);
    %
    ys = round(yn + h * round(A(5,1) * k1 + round(A(5,2) * k2 + round(A(5,3) * k3 + A(5,4) * k4, eps_tt), eps_tt), eps_tt), eps_tt);
    k5 = odefun(tn + c(5) * h, ys);
    %
    ys = yn + h * round(A(6,1) * k1 + round(A(6,2) * k2 + round(A(6,3) * k3 + round(A(6,4) * k4 + A(6,5) * k5, eps_tt), eps_tt), eps_tt), eps_tt);
    k6 = odefun(tn + c(6) * h, ys);
    %
    dy = round(b4(1) * k1 + round(b4(2) * k2 + round(b4(3) * k3 + round(b4(4) * k4 + round(b4(5) * k5 + b4(6) * k6, eps_tt), eps_tt), eps_tt), eps_tt), eps_tt);
    y4 = round(yn + h*dy, eps_tt);
    %
    dy = round(b5(1) * k1 + round(b5(2) * k2 + round(b5(3) * k3 + round(b5(4) * k4 + round(b5(5) * k5 + round(b5(6) * k6 + b5(7) * odefun(tn + h, y4), eps_tt), eps_tt), eps_tt), eps_tt), eps_tt), eps_tt);
    y5 = round(yn + h*dy, eps_tt);
    %
    error_estimate = norm(y5 - y4)/norm(yn);
    %
    if (error_estimate < tol)
      %
      t = tn + h;
      %
      y = y5; 
      %
    end
    %
    %fprintf("t(%d)=%e\n",tstep,t);
    %
    %fprintf("t(%d)=%e\tmax(rank)=%d\n",tstep,t,max(y.r));

    h = h * 0.9 * (tol / error_estimate)^(1/5);  
    %
  end
  %
end
