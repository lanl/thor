function data = explicitSolve(data)
    %
    data = calculateTimeStep(data); 
    %
    rank = zeros(2+data.tt.Neq,10000);
    %
    % save initial ranks
    %
    rank(1,1) = data.curr_time;
    rank(2,1) = data.tt.eps;
    % 
    for eq=1:data.tt.Neq
        rank(eq+2,1) = data.tt.Q{eq}.r(2);
    end
    %
    % advance in time until final time
    %
    if(data.isConsErr)
        %
        Q = remove_ghost(data.tt.Q,data.tt.Neq);
        %
        data.consInt = zeros(3,10000);
        %
        data.consInt(2,1) = sum(Q{1});
        %
        data.consInt(3,1) = sum(fun_KE(Q,data.tt.eps));
        %
    end
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
        data = calculateTimeStep(data); 
        %
        % save ranks
        %
        rank(1,loc_time_step) = data.curr_time;
        rank(2,loc_time_step) = data.tt.eps;
        % 
        for eq=1:data.tt.Neq
            rank(eq+2,loc_time_step) = data.tt.Q{eq}.r(2);
        end
        %
        % save conserved variable integrals
        %
        if(data.isConsErr)
            %
            Q = remove_ghost(data.tt.Q,data.tt.Neq);
            %
            data.consInt(1,loc_time_step) = data.curr_time*data.ref.t;
            %
            data.consInt(2,loc_time_step) =  sum(Q{1});
            %
            data.consInt(3,loc_time_step) = sum(fun_KE(Q,data.tt.eps));
            %
        end
        %
        loc_time_step = loc_time_step + 1;
        %
    end
    %
    data.tt.rank = rank(:,1:loc_time_step-1);
    %
    if(data.isConsErr)
        data.consInt = data.consInt(:,1:loc_time_step-1)';
    else
        data.consInt = [];
    end
    %
    print2screen(data,0,true);
    %
end
%
function KE = fun_KE(Q,eps_tt)
    %
    h  = Q{1};
    hU = Q{2};
    hV = Q{3};
    %
    invh = taylor_inverse(h,eps_tt,100);
    %
    U2 = round(hU.*invh,eps_tt);
    V2 = round(hV.*invh,eps_tt);
    %
    U2 = round(U2.^2,eps_tt);
    V2 = round(V2.^2,eps_tt);
    %
    KE = 0.5*round(U2 + V2,eps_tt);
    %
end
