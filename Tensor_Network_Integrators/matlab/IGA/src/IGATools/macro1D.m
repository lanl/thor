function m1D = macro1D(knot)
    k = ReduceKnotFull(knot);
    len1 = length(k) - 1;
    m1D = zeros(len1, 2);
    % Loop through and populate the cell array
    for i = 1:len1
        m1D(i, :) = [k(i) k(i+1)];
    end
end
