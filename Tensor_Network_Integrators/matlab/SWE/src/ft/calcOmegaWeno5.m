function [Omega0,Omega1,Omega2] = calcOmegaWeno5(params,c0,c1,c2)
    %
    Omega0 = 1./(1.0 + (c1/c0)*params{1,1} + (c2/c0)*params{1,2});  % 1.0 + (c1/c0)*t11 + (c2/c0)*t12;
    Omega1 = 1./(1.0 + (c0/c1)*params{2,1} + (c2/c1)*params{2,2});  % 1.0 + (c0/c1)*t21 + (c2/c1)*t22;
    Omega2 = 1./(1.0 + (c0/c2)*params{3,1} + (c1/c2)*params{3,2});  % 1.0 + (c0/c2)*t31 + (c1/c2)*t32;