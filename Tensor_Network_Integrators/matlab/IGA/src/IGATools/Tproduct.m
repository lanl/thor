function c = Tproduct(a,b)
c = zeros(length(a),length(b));
for i = 1:length(a)
    for j = 1:length(b)
        c(i,j) = a(i)*b(j);
    end 
end
