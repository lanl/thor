function NN = GetNUBSf(knot,zita)
KnotSpan=GetSpan(knot,zita);
[n,p]=EvaluateKnot(knot);
NN   =zeros(1,p+1);
    
    left  = zeros(1,p+1);
    right = zeros(1,p+1);
    
    N_i    = zeros(1,p+1);
    N_i(1) = 1;
    for j= 1 : p
        left(j+1)  = zita-knot(KnotSpan+1-j);
        right(j+1) = knot(KnotSpan+j)-zita;
        saved = 0;
        for r = 0:j-1
            temp     = N_i(r+1)/(right(r+2)+left(j-r+1));
            N_i(r+1) = saved+right(r+2)*temp;
            saved    = left(j-r+1)*temp;
        end
        N_i(j+1) = saved;
    
    end
    NN(:) = N_i;
end