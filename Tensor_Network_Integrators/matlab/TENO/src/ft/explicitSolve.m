function data=explicitSolve(data)
        
    data=calculateTimeStep(data);
    %
    if(data.isConsErr)
        data.consErr = zeros(data.Neq,10000);
    end
    % advance in time until final time
    while(data.curr_time<data.Tfinal)
        %
        fprintf("Time Step= %d\t Time=%.5f\t dt=%.5e\t max_eig=%.2e\n",data.time_step,data.curr_time,data.dt,data.max_eig);
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
        data=calculateTimeStep(data);
        %
    end
    %
    if(data.isConsErr)
        data.consErr = data.consErr(:,1:data.time_step)';
    end
    %
    fprintf("Time Step= %d\t Time=%.5f\t dt=%.5e\t max_eig=%.2e\n",data.time_step,data.curr_time,data.dt,data.max_eig);
    %
end
