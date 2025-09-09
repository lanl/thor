function F=F_nonlinearSWE(Q,data,d)
    %
    g = data.gacc;
    %
    F = cell(1,3);
    %
    h  = Q{1};
    hU = Q{2};
    hV = Q{3};
    % 
    invh = taylor_inverse(h,data.tt.eps,100);
    %
    U = round(hU.*invh,   data.tt.eps);
    V = round(hV.*invh,   data.tt.eps);
    p = round(0.5*g*h.^2, data.tt.eps);
    %
    if(d==1)
        %
        F{1} = hU;
        F{2} = round(hU.*U,    data.tt.eps);
        F{2} = round(F{2} + p, data.tt.eps);
        F{3} = round(hV.*U,    data.tt.eps);
        %
    elseif(d==2)
        %
        F{1} = hV;
        F{2} = round(hU.*V,    data.tt.eps);
        F{3} = round(hV.*V,    data.tt.eps);
        F{3} = round(F{3} + p, data.tt.eps);
        %
    end
    %
end