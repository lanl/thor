function data=forwardEuler(data)
    % get size
    Nx_total = data.Nx_total;
    Ny_total = data.Ny_total;
    Nz_total = data.Nz_total;

    % setup interior indices
    i = [4:Nx_total-3];
    j = [4:Ny_total-3];
    k = [4:Nz_total-3];

    
    % initialize the right-hand side
    data.RHS(:) = 0.0;

    % loop dimensions to set up the right-hand side vector
    for d=1:3
        % calculate shifted indices
        im = i - (1==d);
        jm = j - (2==d);
        km = k - (3==d);
        %
        ip = i + (1==d);
        jp = j + (2==d);
        kp = k + (3==d);

        % calculate fluxes in the d-direction 
        data = faceFlux(data,d);
            
        % accumulate the right-hand side for each equation
        data.RHS(:,:,:,:) = data.RHS(:,:,:,:) ...
                          - (data.FL(:,ip,jp,kp)- data.FL(:,i,j,k)   )/data.h(d) ...
                          - (data.FR(:,i,j,k)   - data.FR(:,im,jm,km))/data.h(d);
    end
    
    % apply source terms if needed
    if (data.isSource)
        data = sourceTerm(data);
    end 
    
    % calculate the solution at the new time level
    data.Q(:,i,j,k) = data.Q(:,i,j,k) + data.dt * data.RHS(:,:,:,:);
    %
end