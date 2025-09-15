function Qstencil = getWeno5Stencil(Q,d)
    %
    Qstencil = cell(1,5);
    %
    Qstencil{1} = applyShift(Q,d,-2);
    Qstencil{2} = applyShift(Q,d,-1);
    Qstencil{3} = Q;
    Qstencil{4} = applyShift(Q,d,+1);
    Qstencil{5} = applyShift(Q,d,+2);
    %