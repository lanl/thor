function count = Findk(knot,t)

  count = 0;
  for i=1:(length(knot) - 1)
   if (t >= knot(i))&& (t < knot(i + 1)), 
       count = i; 
   end
  end

end

