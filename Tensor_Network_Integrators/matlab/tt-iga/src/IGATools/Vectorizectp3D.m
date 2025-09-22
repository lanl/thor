function [ctpxv, ctpyv, ctpzv] = Vectorizectp3D(ctpxn, ctpyn, ctpzn)
[m1, m2, m3] = size(ctpxn);

ctpxv = zeros(m1*m2*m3,1);
ctpyv = zeros(m1*m2*m3,1);
ctpzv = zeros(m1*m2*m3,1);

for i1 = 1:m1
    for i2 = 1:m2
        for i3 = 1:m3
            idx = m1*(i2-1) + i1 + m1*m2*(i3-1);
            ctpxv(idx) = ctpxn(i1,i2,i3);
            ctpyv(idx) = ctpyn(i1,i2,i3);
            ctpzv(idx) = ctpzn(i1,i2,i3);
        end
    end
end