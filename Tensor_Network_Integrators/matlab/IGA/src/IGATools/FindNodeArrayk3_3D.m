function [IndexRe, ctpxRe, ctpyRe, ctpzRe] = FindNodeArrayk3_3D(a, ctpxn, ctpyn, ctpzn)
    [n1, n2, ~] = size(ctpxn);

    a_min = min(a);
    a_max = max(a);

    totalNodes = (a_max - a_min + 1) * n1 * n2;
    ctpxRe = zeros(1, totalNodes);
    ctpyRe = zeros(1, totalNodes);
    ctpzRe = zeros(1, totalNodes);
    IndexRe = zeros(1, totalNodes);
    for i1 = 1:n1
        for i2 = 1:n2
            for i3 = a_min:a_max
                idx = (i1 - 1) + n1 * (i2 - 1) + n1 * n2 * (i3 - a_min) + 1;

                ctpxRe(idx) = ctpxn(i1, i2, i3);
                ctpyRe(idx) = ctpyn(i1, i2, i3);
                ctpzRe(idx) = ctpzn(i1, i2, i3);

                IndexRe(idx) = i1 + n1*(i2 - 1) + n1*n2*(i3 - 1);
            end
        end
    end
end
