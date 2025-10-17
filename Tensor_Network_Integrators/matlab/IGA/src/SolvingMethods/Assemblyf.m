function ft = Assemblyf(fe, InBC)
m = max(InBC(:));
[n1, n2] = size(InBC);
ft = zeros(m,1);
for i=1:n1
    for j=1:n2
        if InBC(i,j)~=-1
            ft(InBC(i,j)) = ft(InBC(i,j)) + fe(j,i);
        end
    end
end
end