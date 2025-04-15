function data=forwardEuler(data)
    % get size
    Nx_total = data.Nx_total;
    Ny_total = data.Ny_total;

    % setup index arrays
    i = [4:Nx_total-3];
    j = [4:Ny_total-3];

    % initialize the right-hand side
    data.RHS(:) = 0.0;

    % loop dimensions to set up the right-hand side vector
    for d=1:2
        % calculate shift indices
        im = i - (1==d);
        jm = j - (2==d);
        %
        ip = i + (1==d);
        jp = j + (2==d); 

        % calculate fluxes in the d-direction 
        data = faceFlux(data,d);

        % accumulate the right-hand side for each equation
        data.RHS(:,:,:) = data.RHS(:,:,:) ...
                        - data.area(d)*(data.FL(:,ip,jp) - data.FL(:,i,j)) ...
                        - data.area(d)*(data.FR(:,i,j)   - data.FR(:,im,jm));
    end
    %
    % add coriolis source term
    %
    data.RHS(2,:,:) = data.RHS(2,:,:) + (data.vol*data.f)*data.Q(3,i,j);
    %
    data.RHS(3,:,:) = data.RHS(3,:,:) - (data.vol*data.f)*data.Q(2,i,j);
    
    % add source terms due to manufactured solutions if needed
    if (data.isManuf)
        data.source(:,:,:) = applyQuadrature(@manuf,data,data.rk_time); % source term is integrated over the cell volumes
        data.RHS(:)        = data.RHS(:) + data.source(:);
    end 
    
    % calculate the solution at the new time level
    data.Q(:,i,j) = data.Q(:,i,j) + (data.dt / data.vol) * data.RHS(:,:,:);
    %
end
