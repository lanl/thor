function Knot = ReduceKnotFull(knot)
  l = length(knot);
  count = 1;
  KnotReduce = zeros(1,l);
  KnotReduce(1) = knot(1);
  for i= 1:(l - 1)
   
  if  knot(i + 1) ~= knot(i)
    
      count = count + 1;
      KnotReduce(count) = knot(i + 1);
  end
  end
  Knot = KnotReduce(1:count);

end

