function m2D = Getmacro2D(knot1,knot2)
  k1 = ReduceKnotFull(knot1);
  k2 = ReduceKnotFull(knot2);
  
  m1 = length(k1);
  m2 = length(k2);
  
  m2D = zeros((m1-1)*(m2-1), 4);
  for i = 1:(m1-1)
   for j = 1:(m2-1)
     m2D((m1-1)*(j-1) + i,:) = [k1(i), k1(i + 1), k2(j), k2(j + 1)];
   end
  end
end

