function data = Q0(data)
    %
    data.Q(:,4:end-3,4:end-3) = applyQuadrature(@solExact,data,0)/data.vol;
    %
    data = data.BC(data, 0);
    %
end