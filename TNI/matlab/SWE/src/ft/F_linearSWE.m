function F = F_linearSWE(Q,dir,g,H)
    %
    % allocate memory
    F = zeros(size(Q));
    %
    eta  = Q(1,:);
    u    = Q(2,:);
    v    = Q(3,:);
    %
    if(dir==1)
        %
        F(1,:) = H*u;
        F(2,:) = g*eta;
        F(3,:) = 0;
        %
    elseif(dir==2)
        %
        F(1,:) = H*v;
        F(2,:) = 0;
        F(3,:) = g*eta;
        %
    end
    %
end
