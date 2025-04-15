function Qrecon=fun_weno5GQ(Qj,weno,nq)
    %--------------------------------------%
    % side=-1 is left
    % side=+1 is right
    %--------------------------------------%
    % 
    weps = weno.weps;
    ngam = size(weno.gamma,2);
    %
    % Get solution in this stencil  
    %
    Qmm = Qj(:,1);
    Qm  = Qj(:,2);
    Q   = Qj(:,3);
    Qp  = Qj(:,4);
    Qpp = Qj(:,5);
    %
    t11 = (  Q   - 2*Qp + Qpp).^2;
    t12 = (3*Q   - 4*Qp + Qpp).^2;
    t21 = (  Qm  - 2*Q  + Qp ).^2;
    t22 = (Qm - Qp).^2         ;
    t31 = (Qmm - 2*Qm + Q    ).^2;
    t32 = (Qmm - 4*Qm + 3*Q  ).^2 ;
    %
    d_13_ov_12 = 13/12;
    d_1_ov_4 = 0.25;
    %
    % compute smoothness indicator
    %
    beta0 = (d_13_ov_12 * t11 + d_1_ov_4 * t12 + weps).^2;
    beta1 = (d_13_ov_12 * t21 + d_1_ov_4 * t22 + weps).^2;
    beta2 = (d_13_ov_12 * t31 + d_1_ov_4 * t32 + weps).^2;
    %
    beta12 = beta1.*beta2;
    beta02 = beta0.*beta2;
    beta01 = beta0.*beta1;
    %
    Qrecon = zeros(size(Q,1),nq);
    %
    
    for ii=1:ngam
        %
        c0 = weno.gamma(1,ii); 
        c1 = weno.gamma(2,ii); 
        c2 = weno.gamma(3,ii);
        %
        c = weno.C{ii};
        %
        Q0 = c(1,1)*Q   + c(1,2)*Qp +  c(1,3)*Qpp;
        Q1 = c(2,1)*Qm  + c(2,2)*Q  +  c(2,3)*Qp ;
        Q2 = c(3,1)*Qmm + c(3,2)*Qm +  c(3,3)*Q  ;
        %
        % calculate weno weights
        %
        om0 = c0*beta12;
        om1 = c1*beta02;
        om2 = c2*beta01;
        %
        denom = om0 + om1 + om2;
        %
        % calculate the reconstructed variable
        %
        if(ii<=nq)
            %
            Qrecon(:,ii) = (Q0.*om0 + Q1.*om1 + Q2.*om2)./denom;
            %
        else
            %
            Qq = (Q0.*om0 + Q1.*om1 + Q2.*om2)./denom;
            %
            Qrecon(:,weno.imid) = weno.sigma_p*Qrecon(:,weno.imid) - weno.sigma_m*Qq;
            %
        end
    end
    %
    return;  
end