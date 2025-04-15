function data = sourceTerm(data)
    %
    % add source terms to the RHS following 
    % Ref: Anti-diffusive flux corrections for high order finite difference WENO schemes
    %      Xu and Shu (2005), Journal of Computational Physics, 
    %      https://doi.org/10.1016/j.jcp.2004.11.014
    %
   data.tt.RHS{3} = round(data.tt.RHS{3} + data.tt.Q{1},data.tt.eps);
   data.tt.RHS{5} = round(data.tt.RHS{5} + data.tt.Q{3},data.tt.eps);
   %
end 