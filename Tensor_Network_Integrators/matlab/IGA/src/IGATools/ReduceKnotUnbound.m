function kub = ReduceKnotUnbound(knot)
k=ReduceKnotFull(knot); 
kub = k(2:(length(k)-1));
end

