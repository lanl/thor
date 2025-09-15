function data = explicitSolve(data)
    %
    data = calculateTimeStep(data);
    %
    % advance in time until final time
    %
    if(data.isConsErr)
        %
        consInt = zeros(3,10000);
        %
        Q=data.Q(:,4:end-3,4:end-3);
        %
        consInt(2,1) = sum(Q(1,:));
        %
        consInt(3,1) = fun_KE(Q);
        %
    end
    %
    loc_time_step = 2;
    %
    while(data.curr_time<data.Tfinal)
        %
        fprintf("Time Step= %d\t Time=%.5f\t dt=%.5e\t max_eig=%.2e\n",data.time_step,data.ref.t*data.curr_time,data.ref.t*data.dt,data.ref.U*data.max_eig);
        %
        % perform rk3
        %
        data = RK3(data);
        %
        % set current time and time step
        %
        data.curr_time = data.curr_time + data.dt;
        data.time_step = data.time_step + 1;
        %
        data = calculateTimeStep(data);
        %
        if(data.isConsErr)
            %
            consInt(1,loc_time_step) = data.curr_time*data.ref.t;
            %
            Q = data.Q(:,4:end-3,4:end-3);
            %
            consInt(2,loc_time_step) = sum(Q(1,:));
            %
            consInt(3,loc_time_step) = fun_KE(Q);
            %
        end
        %
        loc_time_step = loc_time_step + 1;
        %
    end
    %
    fprintf("Time Step= %d\t Time=%.5f\t dt=%.5e\t max_eig=%.2e\n",data.time_step,data.ref.t*data.curr_time,data.ref.t*data.dt,data.ref.U*data.max_eig);
    %
    if(data.isConsErr)
        data.consInt = consInt(:,1:loc_time_step-1)';
    else
        data.consInt = [];
    end
    %
end
%
function KE = fun_KE(Q)
    %
    h  = Q(1,:);
    hU = Q(2,:);
    hV = Q(3,:);
    %
    U2 = hU./h;
    V2 = hV./h;
    %
    U2 = U2.^2;
    V2 = V2.^2;
    %
    KE = sum(0.5*(U2 + V2));
    %
end
