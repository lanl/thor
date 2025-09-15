function params = paramsWeno5(beta0,beta1,beta2);
    %
    params = cell(3,2);
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
    %
end