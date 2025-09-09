function inv_tt = tt_inv(tt,eps_tt,eps_cr)
    %
    dof = sqrt(prod(tt.n));
    %
    nrm = norm(tt)/dof;
    %
    inv_tt  = 2 - tt/nrm;
    %
    inv_tt = round(inv_tt,eps_tt);
    %
    inv_tt = inv_tt/nrm;
    %
    inv_tt  = cross_interpolation({tt}, @(x) fun_inv(x), eps_cr, inv_tt, 0);
    %
end  
%
function y = fun_inv(x)
    y = 1./x(:,1);
end
%
  