function [wnx,Pxn,Pyn,knotyn] = KnotInsertMWY(w,Px,Py,knoty,t)
 n   = length(t);
 wnx = w;
 Pxn = Px;
 Pyn = Py;
 knotyn = knoty;
 for l = 1:n
 [wnx,Pxn,Pyn,knotyn] = KnotInsertWY(wnx,Pxn,Pyn,knotyn,t(l));
 end
 function [wn,Pxn,Pyn,knotyn] = KnotInsertWY(w,Px,Py,knoty,t)
Pwx = Px.*w;
Pwy = Py.*w;
k      = Findk(knoty,t);
[m, p] = EvaluateKnot(knoty);
[n1,n2]= size(w);
wn     = zeros(n1,n2+1);
Pxn    = zeros(n1,n2+1);
Pyn    = zeros(n1,n2+1);
for i = 1:(k - p)
 wn(:,i)    = w(:,i);
 Pxn(:,i)   = Pwx(:,i)./wn(:,i); 
 Pyn(:,i)   = Pwy(:,i)./wn(:,i); 
end

for i = (k - p + 1):k
 alpha     = (t - knoty(i))/(knoty(i + p) - knoty(i));
 wn(:,i)   = (1-alpha)*w(:,i-1)+alpha*w(:,i);
 Pxn(:,i)  = (alpha*Pwx(:,i)+(1-alpha)*Pwx(:,i-1))./wn(:,i);
 Pyn(:,i)  = (alpha*Pwy(:,i)+(1-alpha)*Pwy(:,i-1))./wn(:,i);
end

for i = (k + 1):(m + 1)
 wn(:,i) = w(:,i-1);
 Pxn(:,i) = Pwx(:,i-1)./w(:,i-1);
 Pyn(:,i) = Pwy(:,i-1)./w(:,i-1);
end

knotyn=sort([knoty,t]);
end

end