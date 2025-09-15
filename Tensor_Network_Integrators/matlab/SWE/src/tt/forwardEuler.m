function data=forwardEuler(data)
    % initialize the right-hand side
    for i=1:data.tt.Neq
        data.tt.RHS{i} = data.tt.zero;
    end
    % For upwind schemes, reconstruct quad points over the cells
    if(data.ReconType=="Upwind3")
        data = recon3GQideal(data);
    elseif(data.ReconType=="Upwind5")
        data = recon5GQideal(data);
    end

    %loop dimensions to set up the right-hand side vector
    for d=1:2
        % calculate fluxes in the d-direction 
        data = faceFlux(data,d);
        
        % loop over all equations 
        for i=1:data.tt.Neq
            % accumulate the right-hand side for each equation
            data.tt.RHS{i} = data.tt.RHS{i} ...
                           + data.area(d)*shiftDiff(data.tt.FL{i},d,+1) ...
                           - data.area(d)*shiftDiff(data.tt.FR{i},d,-1) ;
        end
        %
    end
    %
    % add coriolis source term
    %
    data.tt.RHS{2} = data.tt.RHS{2} + (data.vol*data.f)*data.tt.Q{3};
    %
    data.tt.RHS{3} = data.tt.RHS{3} - (data.vol*data.f)*data.tt.Q{2};
    
    % add source terms due to manufactured solutions if needed
    if (data.isManuf)
        %
        data = manuf(data.tt.q_xy{1},data.tt.q_xy{2},data.rk_time,data);
        %
        for i=1:data.tt.Neq
            %
            data.tt.RHS{i} = round(data.tt.RHS{i} + data.tt.source{i},data.tt.eps_rk(i));
            %
        end
        %
    end 
    %
    for d=1:2
        data.tt.RHS = zero_ghost(data.tt.RHS,d);
    end 
    
    % calculate the solution at the new time level
    for i=1:data.tt.Neq
        %
        data.tt.Q{i} = data.tt.Q{i} + (data.dt / data.vol) * data.tt.RHS{i};
        data.tt.Q{i} = round(data.tt.Q{i},data.tt.eps_rk(i));
        %
    end
    %
end