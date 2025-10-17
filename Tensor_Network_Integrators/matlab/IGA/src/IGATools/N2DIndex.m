function [INX, m2D] = N2DIndex(knot1,knot2)
N2DDomain = GetNURBS2DDomain(knot1,knot2);
m2D       = Getmacro2D(knot1,knot2);
[NOE, ~]  = size(m2D);
[L, ~]    = size(N2DDomain);
[~, p1]   = EvaluateKnot(knot1);
[~, p2]   = EvaluateKnot(knot2);
INX       = zeros(NOE,(p1+1)*(p2+1));
disp('Compute connectivity');
tic
for k=1:NOE
   count=0;
   for i=1:L
    if ~(N2DDomain(i,2) <=  m2D(k,1)||...
         N2DDomain(i,1) >=  m2D(k,2)||...
         N2DDomain(i,4) <=  m2D(k,3)||...    
         N2DDomain(i,3) >=  m2D(k,4))
     count=count+1;
     INX(k,count) = i;
    end
    if count==(p1+1)*(p2+1)
        break
    end
   end
end 
toc  
end