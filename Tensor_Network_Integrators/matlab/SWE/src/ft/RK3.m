function data = RK3(data)
    % save solution at the previous time step
    Qprev = data.Q;
    %
    % loop over rk stages
    %
    for rki=1:data.rkmax
        data.rki = rki;
        % get rk coeffcients
        c1 = data.crk(rki,1);
        c2 = data.crk(rki,2);
        % get rk time step coefficient
        cdt = data.cdt(rki);
        % calculate rk time
        data.rk_time = data.curr_time + cdt*data.dt;
        % forward euler step
        data = forwardEuler(data);
        % update solution
        data.Q(:) = c1*Qprev(:) + c2*data.Q(:);
        % apply boundary conditions         
        if(rki<data.rkmax)
            bc_time = data.curr_time + data.cdt(rki+1)*data.dt;
        else
            bc_time = data.curr_time + data.dt;
        end
        data = data.BC(data, bc_time);
        %
    end
    %
end
