function NURBS = GetNURBSf(knot,zita,w)
NUBS  = GetNUBSf(knot,zita);
NURBS = NUBS/sum(w.*NUBS);
end

