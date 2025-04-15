function data = RK3(data)
    % save solution at the previous time step
    Qprev = data.Q;
    %
    if(data.isConsErr)
        %
        time_step = data.time_step + 1;
        %
        RHS = zeros(size(data.RHS));
        %
        coeff = [1/6;1/6;2/3];
        %
    end
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
        data=forwardEuler(data);
        % update solution
        data.Q(:) = c1*Qprev(:) + c2*data.Q(:);
        % apply boundary conditions         
        if(rki<data.rkmax)
            bc_time = data.curr_time + data.cdt(rki+1)*data.dt;
        else
            bc_time = data.curr_time + data.dt;
        end
        data.Q(:,:,:,:) = data.BC(data.Q,data.X,data.Y,data.Z,bc_time,data.gam);

        if(data.isConsErr)
            RHS(:) = RHS(:) + (coeff(rki) * data.dt) * data.RHS(:);
        end
        %
    end
    %
    if(data.isConsErr)
        i = [4:data.Nx_total-3];
        j = [4:data.Ny_total-3];
        k = [4:data.Nz_total-3];
        %
        dQ = data.Q(:,i,j,k) - Qprev(:,i,j,k) - RHS;
        %
        for eq=1:data.Neq
            %
            data.consErr(eq,time_step) = max(abs(dQ(eq,:,:,:)),[],"all");
            %
        end
    end
    %
end
