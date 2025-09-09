function [y,yerr] = ERKstep_embedded(fcn, y0, t0, dt, B)
% Inputs:
%    fcn = ODE RHS function, f(t,y)
%    y0  = solution at beginning of time step
%    t0  = 'time' at beginning of time step
%    dt   = step size to take
%    B   = Butcher table to use
%
% Outputs:
%    y     = new solution at t0+dt
%    yerr  = error vector

% extract ERK method information from B
[Brows, Bcols] = size(B);
s = Bcols - 1;        % number of stages
c = B(1:s,1);         % stage time fraction array
b = (B(s+1,2:s+1))';  % solution weights (convert to column)
A = B(1:s,2:s+1);     % RK coefficients
d = (B(s+2,2:s+1))';  % embedding coefficients

% initialize storage for RHS 
k = cell(1,s);

% loop over stages
for stage=1:s

  % construct stage solution and evaluate RHS
  %    zi = y_n + dt*sum_{j=1}^{i-1} (A(i,j)*f(zj))
  z = y0;
  for j=1:stage-1
     z = z + dt*A(stage,j)*k{j};
  end

  % construct new stage RHS
  % Assume truncation happens in fcn 
  k{stage} = fcn(t0+dt*c(stage),z);

end

% compute new solution and error estimate
%    ynew = yold + dt*sum(b(j)*fj)
% kb = cellfun(@times,k,num2cell(dt*b),'UniformOutput',false);
y = y0;
yerr = 0;
for i=1:s
   y = y + dt*b(i)*k{i};
   yerr = yerr + dt*(b(i)-d(i))*k{i};
end

end % end of function