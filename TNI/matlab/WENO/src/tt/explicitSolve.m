function data=explicitSolve(data)
    %
    data=calculateTimeStep(data); 
    %
    rank = zeros(2+2*data.tt.Neq,10000);
    %
    if(data.isConsErr)
        data.consErr = zeros(data.tt.Neq,10000);
    end
    %
    % save initial ranks
    %
    rank(1,1) = data.curr_time;
    rank(2,1) = data.tt.eps;
    % 
    for eq=1:data.tt.Neq
        rank(2*eq+1:2*eq+2,1) = data.tt.Q{eq}.r(2:3);
    end
    %
    % advance in time until final time
    %
    loc_time_step = 2;
    %
    while(data.curr_time<data.Tfinal)
        % 
        print2screen(data,0,false);
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
        % save ranks
        %
        rank(1,loc_time_step) = data.curr_time;
        rank(2,loc_time_step) = data.tt.eps;
        % 
        for eq=1:data.tt.Neq
            rank(2*eq+1:2*eq+2,loc_time_step) = data.tt.Q{eq}.r(2:3);
        end
        %
        loc_time_step = loc_time_step + 1;
        %
    end
    %
    if(data.isConsErr)
        data.consErr = data.consErr(:,1:data.time_step)';
    else
        data.consErr = [];
    end
    %
    data.tt.rank = rank(:,1:loc_time_step-1);
    %
    print2screen(data,0,true);
    %
end
