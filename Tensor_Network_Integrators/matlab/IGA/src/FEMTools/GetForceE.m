function [fexte, PM] = GetForceE(IndexEBC, BC, PIndex, PValues)
[n1, n2] = size(IndexEBC);
PM = PIndex;
for j = 1:length(PIndex)
    for i = 1:length(BC)
        if BC(i) < PIndex(j)
            PM(j) = PM(j)-1;
        end
    end
end
fexte = zeros(n2,n1);

for k = 1:length(PIndex)
    for i = 1:n1
        for j = 1:n2
            if IndexEBC(i,j) == PM(k)
                fexte(j,i) = PValues(k);
            end
        end
    end
end