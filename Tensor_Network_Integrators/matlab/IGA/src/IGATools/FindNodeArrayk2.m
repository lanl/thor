function [IndexRe, ctpxRe, ctpyRe] = FindNodeArrayk2(a, ctpxn, ctpyn)
    [n1, n2] = size(ctpxn);
    
    a_min = min(a);
    a_max = max(a);
    
    totalElements = n1 * (a_max - a_min + 1);
    ctpxRe = zeros(1, totalElements);
    ctpyRe = zeros(1, totalElements);
    IndexRe = zeros(1, totalElements);
    
    for i = 1:n1
        for j = a_min:a_max
            idx = (j - a_min) * n1 + i;
            ctpxRe(idx) = ctpxn(i, j);
            ctpyRe(idx) = ctpyn(i, j);
            IndexRe(idx) = (j - 1) * n1 + i;
        end
    end
end