function params = paramsWeno5(Qmm,Qm,Q,Qp,Qpp,weps);
    %
    params = cell(3,2);
    %
    % compute smoothness parameter
    %
    beta0 = (13/12*( Q  - 2*Qp + Qpp).^2 + 1/4*(3*Q - 4*Qp + Qpp).^2 + weps);
    %
    beta1 = (13/12*(Qm  -  2*Q +  Qp).^2 + 1/4*(Qm - Qp).^2 + weps);
    %
    beta2 = (13/12*(Qmm - 2*Qm +   Q).^2 + 1/4*(Qmm - 4*Qm + 3*Q).^2 + weps);
    %
    wp = 2;
    %
    params{1,1} = (beta0./beta1).^wp;  % t11
    params{2,1} = (beta1./beta0).^wp;  % t21
    params{3,1} = (beta2./beta0).^wp;  % t31
    %
    params{1,2} = (beta0./beta2).^wp;  % t12
    params{2,2} = (beta1./beta2).^wp;  % t22
    params{3,2} = (beta2./beta1).^wp;  % t32