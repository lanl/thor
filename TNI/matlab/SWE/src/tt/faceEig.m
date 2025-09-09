function data=faceEig(data,d)
    % 
    if(data.SWEtype=="linear")
        %
        Q = data.tt.QL; % eig is just a constant value for the linear equation. So, Q value is not important
        %
    else
        %
        % calculate the average value at the face
        %
        Q    = cell(1,data.tt.Neq);
        %
        for i = 1:data.tt.Neq
            %
            Q{i} = 0.5*round(data.tt.QL{i} + applyShift(data.tt.QR{i},d,-1),data.tt.eps_rk(i));
            %
        end
        %
    end
    %
    data.tt.eig = data.Eig(Q,data,d);
    %
end