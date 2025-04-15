function data = RK3(data)
    % save solution at the previous time step
    Qprev = data.tt.Q;
    %
    % loop over rk stages
    %
    for rki=1:data.rkmax
        %
        print2screen(data,rki,false);
        
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
        for i=1:data.tt.Neq
            data.tt.Q{i} = c2*data.tt.Q{i} + c1*Qprev{i};
            data.tt.Q{i} = round(data.tt.Q{i}, data.tt.eps_rk(i));
        end
        % apply boundary conditions      
        if(rki<data.rkmax)
            bc_time = data.curr_time + data.cdt(rki+1)*data.dt;
        else
            bc_time = data.curr_time + data.dt;
        end
        data.tt.Q = applyBC(data.tt.Q,data,bc_time);
        %
    end
end