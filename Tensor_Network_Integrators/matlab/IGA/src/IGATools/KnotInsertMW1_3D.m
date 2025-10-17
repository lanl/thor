function [wn, Pxn, Pyn, Pzn, knotxn] = KnotInsertMW1_3D(w, Px, Py, Pz, knotx, t)
 n   = length(t);
 wn = w;
 Pxn = Px;
 Pyn = Py;
 Pzn = Pz;
 knotxn = knotx;
 for l = 1:n
 [wn, Pxn, Pyn, Pzn, knotxn] = KnotInsertW1(wn, Pxn, Pyn, Pzn, knotxn, t(l));
 end
 %% Function for multiple time refinement on 1st-direction
 function [wn, Pxn, Pyn, Pzn, knotxn] = KnotInsertW1(w, Px, Py, Pz, knotx, t)
     Pwx = Px.*w;
     Pwy = Py.*w;
     Pwz = Pz.*w;
     k      = Findk(knotx, t);
     [m, p] = EvaluateKnot(knotx);
     [n1, n2, n3] = size(w);
     wn     = zeros(n1+1, n2, n3);
     Pxn    = zeros(n1+1, n2, n3);
     Pyn    = zeros(n1+1, n2, n3);
     Pzn    = zeros(n1+1, n2, n3);
     for i = 1:(k - p)
         wn(i,:,:)    = w(i,:,:);
         Pxn(i,:,:)   = Pwx(i,:,:)./wn(i,:,:); 
         Pyn(i,:,:)   = Pwy(i,:,:)./wn(i,:,:); 
         Pzn(i,:,:)   = Pwz(i,:,:)./wn(i,:,:); 
     end
     for i = (k - p + 1):k
         alpha      = (t - knotx(i))/(knotx(i + p) - knotx(i));
         wn(i,:,:)  = (1-alpha)*w(i-1,:,:)+alpha*w(i,:,:);
         Pxn(i,:,:) = (alpha*Pwx(i,:,:)+(1-alpha)*Pwx(i-1,:,:))./wn(i,:,:);
         Pyn(i,:,:) = (alpha*Pwy(i,:,:)+(1-alpha)*Pwy(i-1,:,:))./wn(i,:,:);
         Pzn(i,:,:) = (alpha*Pwz(i,:,:)+(1-alpha)*Pwz(i-1,:,:))./wn(i,:,:);
     end
     for i = (k + 1):(m + 1)
         wn(i,:,:)  = w(i-1,:,:);
         Pxn(i,:,:) = Pwx(i-1,:,:)./wn(i,:,:);
         Pyn(i,:,:) = Pwy(i-1,:,:)./wn(i,:,:);
         Pzn(i,:,:) = Pwz(i-1,:,:)./wn(i,:,:);
     end

knotxn = sort([knotx,t]);
end
%%
end