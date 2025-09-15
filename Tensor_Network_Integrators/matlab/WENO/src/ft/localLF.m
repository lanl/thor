function data=localLF(data,d)
    % 
    F = data.F(data.Q,d,data.gam);
    %
    data.invFL(:,:,:,:) = 0.5*(F - data.eig.*data.Q);
    data.invFR(:,:,:,:) = 0.5*(F + data.eig.*data.Q);
    %
end
    