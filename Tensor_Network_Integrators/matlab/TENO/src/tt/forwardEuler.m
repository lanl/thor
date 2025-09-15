function data=forwardEuler(data)
    % initialize the right-hand side
    for i=1:data.tt.Neq
        data.tt.RHS{i} = data.tt.zero;
    end

    %loop dimensions to set up the right-hand side vector
    for d=1:3
        % calculate fluxes in the d-direction 
        data = faceFlux(data,d);
        % loop over all equations 
        for i=1:data.tt.Neq
            % accumulate the right-hand side for each equation
            data.tt.RHS{i} = data.tt.RHS{i} ...
                           + shiftDiff(data.tt.FL{i},d,+1)/data.h(d) ...
                           - shiftDiff(data.tt.FR{i},d,-1)/data.h(d) ;
            %
            data.tt.RHS{i} = round(data.tt.RHS{i},data.tt.eps);
        end
        %
    end
    %
    % apply source terms if needed
    %
    if (data.isSource)
        data = sourceTerm(data);
    end 
    %
    for d=1:3
        data.tt.RHS = zero_ghost(data.tt.RHS,d);
    end 
    
    % calculate the solution at the new time level
    for i=1:data.tt.Neq
        if(data.vol*norm(data.tt.RHS{i})/prod(data.tt.RHS{i}.n)>1e-14)
            data.tt.Q{i} = round(data.tt.Q{i} + data.dt * data.tt.RHS{i},data.tt.eps_rk{i});
        end
    end
    %
end
