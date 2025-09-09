function data=Q0(data)
    %
    gam = data.gam;
    %
    x0    = 0.5;    % initial discontinuity
    p_r   = 0.1;    % right BC for pressure 
    u_r   = 0;      % right BC for velocity
    rho_r = 0.125;  % right BC for density
    p_l   = 1;      % left BC for pressure
    u_l   = 0;      % left BC for velocity
    rho_l = 1;      % left BC for density
    %
    for i=1:size(data.X,1)
        %
        if(data.X(i,1,1)<=x0)
            data.Q(1,i,:,:) = rho_l; 
            data.Q(5,i,:,:) = p_l/(gam-1); 
        else
            data.Q(1,i,:,:) = rho_r; 
            data.Q(5,i,:,:) = p_r/(gam-1); 
        end
        %
    end
    %
end