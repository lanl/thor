function data = RK3(data)
    % save solution at the previous time step
    Qprev = data.tt.Q;
    %
    if(data.isConsErr)
        %
        time_step = data.time_step + 1;
        %
        RHS = cell(1,data.tt.Neq);
        %
        coeff = [1/6;1/6;2/3];
        %
        for eq=1:data.tt.Neq
            RHS{eq} = data.tt.zero;
        end
        %
    end
    %
    % loop over rk stages
    %
    for rki=1:data.rkmax
        %
        print2screen(data,rki,false);
        %
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
            data.tt.Q{i} = round(c1*Qprev{i} + c2*data.tt.Q{i}, data.tt.eps);
        end
        % apply boundary conditions      
        if(rki<data.rkmax)
            bc_time = data.curr_time + data.cdt(rki+1)*data.dt;
        else
            bc_time = data.curr_time + data.dt;
        end
        data.tt.Q=applyBC(data.tt.Q,data,bc_time);
        %
        if(data.isConsErr)
            for eq=1:data.tt.Neq
                RHS{eq} = RHS{eq} + (coeff(rki) * data.dt) * data.tt.RHS{eq};
            end
        end
        %
    end
    %
    if(data.isConsErr)
        %
        dQ = cell(1,data.tt.Neq);
        %
        for eq=1:data.tt.Neq
            %
            dQ{eq} = data.tt.Q{eq}-Qprev{eq}-RHS{eq};
            %
        end
        %
        dQ = remove_ghost(dQ,data.tt.Neq);
        %
        for eq=1:data.tt.Neq
            %
            data.consErr(eq,time_step) = tt_max_abs(dQ{eq});
            %
        end
        %
    end
    %
end
%