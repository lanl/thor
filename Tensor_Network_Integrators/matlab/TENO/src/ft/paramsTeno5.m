function params = paramsTeno5(beta0,beta1,beta2);
    %
    params = cell(1,3);
    %
    Ct   = 1e-5;
    wp   = 6;
    %
    tau5=abs(beta2-beta0);
    %
    gam0 = (1 + tau5./beta0).^wp;
    gam1 = (1 + tau5./beta1).^wp;
    gam2 = (1 + tau5./beta2).^wp;
    %
    sumG = gam0 + gam1 + gam2;
    %
    chi0 = gam0./sumG;
    chi1 = gam1./sumG;
    chi2 = gam2./sumG;
    %
    params{1} = zeros(size(chi0)); % delta0
    params{2} = zeros(size(chi1)); % delta1
    params{3} = zeros(size(chi2)); % delta2
    %
    params{1}(chi0>=Ct) = 1;
    params{2}(chi1>=Ct) = 1;
    params{3}(chi2>=Ct) = 1;
    %
end