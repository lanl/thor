function data=calculateTimeStep(data)
    %
    data=calculateEps(data);
    %
    % calculate max eigenvalue
    %
    if(data.dtOption==1)
        %
        [~,idx] = getMaxRank(data.tt.Q);
        %
        Eig = cross_interpolation(data.tt.Q, @(xxx) data.Eig(xxx,4,data.gam), data.tt.eps_cr, data.tt.Q{idx}, 0);
        %
        Eig = q_unfold(Eig,data.tt.eps);
        %
        data.max_eig = tt_max_abs(Eig{1});
        %
    else    
        data.max_eig = 0.0;
    end
    
    % calculate time step
    if(data.dtOption==1)
        data.dt = data.cfl*min(data.h)/data.max_eig;
    elseif(data.dtOption==2)
        data.dt = data.cfl*(min(data.h)^(5/3));
    else
        data.dt = data.cfl*min(data.h);
    end 
    
    % check dt so that data.curr_time will not exceed Tfinal
    if(data.curr_time + data.dt > data.Tfinal)
        data.dt = data.Tfinal - data.curr_time;
    end
    %
end
