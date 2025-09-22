function C = CrossProduct2D(A,B)
m=length(A);
n=length(B);

C=zeros(m,n);
for i=1:m
    for j=1:n
    C(i,j)=A(i)*B(j);
    end
end

end

