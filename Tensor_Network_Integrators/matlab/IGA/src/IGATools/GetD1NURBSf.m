function D1NURBS = GetD1NURBSf(knot,zita,w)
NN    = GetDNUBSf(knot,zita,1);
NUBS  = NN(:,1);
DNUBS = NN(:,2);

SumNUBS  = sum(w.*NUBS);
SumDNUBS = sum(w.*DNUBS);

NURBS   = NUBS/SumNUBS;
DNURBS  = ((DNUBS.*w)*SumNUBS - (NUBS.*w)*SumDNUBS)/(SumNUBS^2);

D1NURBS =[NURBS DNURBS];
end

