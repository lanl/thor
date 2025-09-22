function d = Tproduct3(a, b, c)
d = zeros(length(a), length(b), length(c));
for i = 1:length(a)
    for j = 1:length(b)
        for k = 1:length(c)
            d(i, j, k) = a(i)*b(j)*c(k);
        end 
    end
end
