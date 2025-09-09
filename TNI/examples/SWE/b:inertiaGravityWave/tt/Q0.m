function data = Q0(data)
    % calculate cell-averaged initial condition
    data.tt.Q = solExact(data.tt.q_xy{1},data.tt.q_xy{2},data.curr_time,data);
end