function data=sourceTerm(data)
    %
    % calculate source term
    % 
    data.source(:,:,:,:) = solManuf(data.X(4:end-3,4:end-3,4:end-3),data.Y(4:end-3,4:end-3,4:end-3),data.Z(4:end-3,4:end-3,4:end-3),data.rk_time,data.gam);
    %
    % update right-hand side
    %
    data.RHS(:) = data.RHS(:) + data.source(:);
    %
end