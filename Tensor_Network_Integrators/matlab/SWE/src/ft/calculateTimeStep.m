function data=calculateTimeStep(data)
    % calculate max eigenvalue
    if(data.dtOption==1 || data.dtOption==2)
        data.max_eig = max(data.Eig(data.Q,3,data.gacc,data.H),[],"all");
    else    
        data.max_eig = 0.0;
    end
    
    % calculate time step
    if(data.dtOption==1)
        data.dt = data.cfl*min(data.h)/data.max_eig;
    elseif(data.dtOption==2)
        data.dt = data.cfl*(min(data.h)^(data.hp/data.rkmax));
    else
        data.dt = data.cfl*min(data.h);
    end 
    % check dt so that data.curr_time will not exceed Tfinal
    if(data.curr_time + data.dt > data.Tfinal)
        data.dt = data.Tfinal - data.curr_time;
    end
    %
end