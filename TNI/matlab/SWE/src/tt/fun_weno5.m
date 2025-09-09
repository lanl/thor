function Qrecon=fun_weno5(Qj,side,weps)
    %--------------------------------------%
    % side=-1 is left
    % side=+1 is right
    %--------------------------------------%
    %
    % Get solution in this stencil  
    %
    if(length(size(Qj))==2)
        Qmm = Qj(:,1);
        Qm  = Qj(:,2);
        Q   = Qj(:,3);
        Qp  = Qj(:,4);
        Qpp = Qj(:,5);
    elseif(length(size(Qj))==3)
        Qmm = Qj(:,:,1);
        Qm  = Qj(:,:,2);
        Q   = Qj(:,:,3);
        Qp  = Qj(:,:,4);
        Qpp = Qj(:,:,5);
    end
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
    % reconstruct interpolation polynomials on the right side or left side 
    %
    if(side==1)
        Q0 = ( 2*Q   + 5*Qp -    Qpp)/6;
        Q1 = (  -Qm  + 5*Q  +  2*Qp )/6;
        Q2 = ( 2*Qmm - 7*Qm + 11*Q  )/6;
    elseif (side==-1)
        Q0 = ( 11*Q   - 7*Qp + 2*Qpp)/6;
        Q1 = (  2*Qm  + 5*Q  -   Qp )/6;
        Q2 = (   -Qmm + 5*Qm + 2*Q  )/6;
    else
        error("side value can be either 1 or -1.");
    end
    %
    if(side==1)
        c0 = 0.3; c1 = 0.6; c2 = 0.1;
    elseif(side==-1)
        c0 = 0.1; c1 = 0.6; c2 = 0.3;
    end
    %
    % calculate weno weights
    %
    om0 = c0*beta1.*beta2;
    om1 = c1*beta0.*beta2;
    om2 = c2*beta0.*beta1;
    %
    denom = om0 + om1 + om2;
    %
    % calculate the reconstructed variable
    %
    Qrecon = (Q0.*om0 + Q1.*om1 + Q2.*om2)./denom;
    %
    return;  
end