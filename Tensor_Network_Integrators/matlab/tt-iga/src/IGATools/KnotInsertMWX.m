function [wnx,Pxn,Pyn,knotxn] = KnotInsertMWX(w,Px,Py,knotx,t)
 n   = length(t);
 wnx = w;
 Pxn = Px;
 Pyn = Py;
 knotxn = knotx;
 for l = 1:n
 [wnx,Pxn,Pyn,knotxn] = KnotInsertWX(wnx,Pxn,Pyn,knotxn,t(l));
 end 
 function [wn,Pxn,Pyn,knotxn] = KnotInsertWX(w,Px,Py,knotx,t)
Pwx = Px.*w;
Pwy = Py.*w;
k      = Findk(knotx,t);
[m, p] = EvaluateKnot(knotx);
[n1,n2]= size(w);
wn     = zeros(n1+1,n2);
Pxn    = zeros(n1+1,n2);
Pyn    = zeros(n1+1,n2);
for i = 1:(k - p)
 wn(i,:)    = w(i,:);
 Pxn(i,:)   = Pwx(i,:)./wn(i,:); 
 Pyn(i,:)   = Pwy(i,:)./wn(i,:); 
end

for i = (k - p + 1):k
 alpha     = (t - knotx(i))/(knotx(i + p) - knotx(i));
 wn(i,:)   = (1-alpha)*w(i-1,:)+alpha*w(i,:);
 Pxn(i,:)  = (alpha*Pwx(i,:)+(1-alpha)*Pwx(i-1,:))./wn(i,:);
 Pyn(i,:)  = (alpha*Pwy(i,:)+(1-alpha)*Pwy(i-1,:))./wn(i,:);
end

for i = (k + 1):(m + 1)
 wn(i,:) = w(i-1,:);
 Pxn(i,:) = Pwx(i-1,:)./wn(i,:);
 Pyn(i,:) = Pwy(i-1,:)./wn(i,:);
end

knotxn = sort([knotx,t]);
end

end