function Kuut=AssemblyK(Kuue,In1,In2)
m=max(In1(:));
n=max(In2(:));
Kuut=zeros(m,n);
[n1,n2]=size(In1);
[~ ,n3]=size(In2);
 for i=1:n1
  for j=1:n2
    for k=1:n3
      if In1(i,j) ~= -1 && In2(i, k) ~= -1
        Kuut(In1(i,j),In2(i,k)) = Kuut(In1(i,j),In2(i,k)) + Kuue(j,k,i);
      end
     end
  end
 end
end