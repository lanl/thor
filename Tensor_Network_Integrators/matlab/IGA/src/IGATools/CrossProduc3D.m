function D = CrossProduc3D(A,B,C)
m=length(A);
n=length(B);
l=length(C);

for i=1:m
    for j=1:n
        for k=1:l
            D(i,j,k)=A(i)*B(j)*C(k);
        end
    end
end
end

