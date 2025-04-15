function Qstencil = getWeno5StencilNeq(Q,d)
    %
    Neq = length(Q);
    %
    Qstencil = cell(5*Neq,1);
    %
    for eq=1:Neq
        for j=1:5
            Qstencil{eq + Neq*(j-1)} = applyShift(Q{eq},d,j-3);
        end
    end
    %
end