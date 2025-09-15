function data = sourceTerm(data)
    %
    % add source terms to the RHS following 
    % Ref: Anti-diffusive flux corrections for high order finite difference WENO schemes
    %      Xu and Shu (2005), Journal of Computational Physics, 
    %      https://doi.org/10.1016/j.jcp.2004.11.014
    %
    data.RHS(3,:,:,:) = data.RHS(3,:,:,:) + data.Q(1,4:end-3,4:end-3,4:end-3);
    data.RHS(5,:,:,:) = data.RHS(5,:,:,:) + data.Q(3,4:end-3,4:end-3,4:end-3);
    %
end