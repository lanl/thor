function Index2D=GetIndex2D(Index1D)
[n1,n2]=size(Index1D);
Index2D=zeros(n1,2*n2);
    for i=1:n1
        temp=[2*Index1D(i,:)-1;2*Index1D(i,:)];
        Index2D(i,:)=temp(:);
    end
end