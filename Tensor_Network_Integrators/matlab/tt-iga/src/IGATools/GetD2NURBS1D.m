function D2NURBS = GetD2NURBS1D(knot,zita,w)
NN     = GetDNUBS(knot,zita,2);
NUBS   = NN(:,1);
D1NUBS = NN(:,2);
D2NUBS = NN(:,3);

SumNUBS   = sum(w.*NUBS);
SumD1NUBS = sum(w.*D1NUBS);
SumD2NUBS = sum(w.*D2NUBS);

NURBS    = NUBS/SumNUBS;

TS1      = (D1NUBS.*w)*SumNUBS - (NUBS.*w)*SumD1NUBS;
MS1      = SumNUBS^2;
D1NURBS  = TS1/MS1;

TS2      = (D2NUBS*SumNUBS+D1NUBS*SumD1NUBS-D1NUBS*SumD1NUBS-NUBS*SumD2NUBS)*MS1...
           -(2*SumNUBS*SumD1NUBS)*TS1;
MS2      = MS1^2;
D2NURBS  = TS2/MS2;


D2NURBS =[NURBS D1NURBS D2NURBS];
end

