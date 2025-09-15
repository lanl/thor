function Q=BCisentropicVortex(Q,data,t,d)
    %
    % To avoid repeated exact BC computations, only compute them for d=1
    %
    if(d==1)
        %
        % set ghost cells to zero
        %
        for j=1:3
            Q = zero_ghost(Q,j);
        end
        %
        % calculate exact solution at t 
        %
        Qbc = solExact(data,t,Q); % Q here is used as an initial guess for the cross interpolation
        %
        % To extract the exact ghost cell values, extract the exact interior cell values
        %
        Qinterior = Qbc;
        %
        for j=1:3
            Qinterior = zero_ghost(Qinterior,j);
        end
        %
        % Subtract interior values from the exact solution to obtain the exact ghost cell values
        %
        for i=1:length(Qinterior)
            %
            % This TT contains only ghost cell values filled with exact solution (interior is zero)
            %
            Qbc{i} = round(Qbc{i} - Qinterior{i},data.tt.eps);
            %
            % Now apply the BC to the numerical solution for this variable 
            %
            Q{i} = Q{i} + Qbc{i};
            %
        end
        %
        % Finally apply periodic BC in the z-direction
        %
        Q=BCperiodic(Q,data,t,3);
        %
    end
    %
end