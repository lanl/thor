function data=localLF(data,d)
    %
    data.tt.FR = data.F(data.tt.QR,data,d);
    data.tt.FL = data.F(data.tt.QL,data,d);
    %
    if(data.SWEtype=="linear")
        eigR = data.tt.eig; % since it is a constant
    else
        eigR = applyShift(data.tt.eig,d,+1);
    end
    %
    for i=1:data.tt.Neq
        %
        eQL = data.tt.eig.*data.tt.QL{i};
        eQR = eigR.*data.tt.QR{i};
        %
        % calculate fluxes at face quadrature points
        data.tt.FL{i} = 0.5*(data.tt.FL{i} - eQL);
        data.tt.FR{i} = 0.5*(data.tt.FR{i} + eQR);
        % apply quadrature rule to calculate surface integral of fluxes for high-order schemes
        if(data.hp > 1)
            data.tt.FL{i} = applyQuadrature(data.tt.FL{i},data.quad,data.h,d);
            data.tt.FR{i} = applyQuadrature(data.tt.FR{i},data.quad,data.h,d);
        end
        %
    end
    %
end