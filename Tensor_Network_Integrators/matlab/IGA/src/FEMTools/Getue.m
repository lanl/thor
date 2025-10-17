function ue = Getue(ut,uInBC) 
[m1, m2] = size(uInBC);
ue = zeros(m1,m2);
  for i=1:m1
   for j=1:m2
    if uInBC(i,j) ~= -1
        ue(i,j) = ut(uInBC(i,j));
    end
   end
  end
  
end