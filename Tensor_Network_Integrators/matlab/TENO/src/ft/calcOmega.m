function [Omega0,Omega1,Omega2] = calcOmega(params,c0,c1,c2)
    %
    sumd =  c0*params{1} + c1*params{2} + c2*params{3}; % c0*delta0 + c1*delta1 + c2*delta2;
    %
    Omega0 = c0*params{1}./sumd; % c0*delta0./sumd;
    Omega1 = c1*params{2}./sumd; % c1*delta1./sumd;
    Omega2 = c2*params{3}./sumd; % c2*delta2./sumd;
    %
end