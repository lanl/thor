function [IndexRe, ctpxRe, ctpyRe] = FindNodeArrayk1(a, ctpxn, ctpyn)
    [n1, n2] = size(ctpxn);
    
    a_min = min(a);
    a_max = max(a);
    
    totalElements = n2 * (a_max - a_min + 1);
    ctpxRe = zeros(1, totalElements);
    ctpyRe = zeros(1, totalElements);
    IndexRe = zeros(1, totalElements);

    for i = a_min:a_max
        for j = 1:n2
            idx = (i - a_min) * n2 + j;
            ctpxRe(idx) = ctpxn(i, j);
            ctpyRe(idx) = ctpyn(i, j);
            IndexRe(idx) = (j - 1) * n1 + i;
        end
    end
end