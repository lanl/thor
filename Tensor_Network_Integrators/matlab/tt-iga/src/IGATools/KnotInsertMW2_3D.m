function [wn, Pxn, Pyn, Pzn, knotyn] = KnotInsertMW2_3D(w, Px, Py, Pz, knoty, t)
 n   = length(t);
 wn = w;
 Pxn = Px;
 Pyn = Py;
 Pzn = Pz;
 knotyn = knoty;
 for l = 1:n
     [wn, Pxn, Pyn, Pzn, knotyn] = KnotInsertW2(wn, Pxn, Pyn, Pzn, knotyn, t(l));
 end
 %% Function for multiple time refinement on 2nd-direction
    function [wn, Pxn, Pyn, Pzn, knotyn] = KnotInsertW2(w, Px, Py, Pz, knoty, t)
        Pwx = Px.*w;
        Pwy = Py.*w;
        Pwz = Pz.*w;
        k      = Findk(knoty, t);
        [m, p] = EvaluateKnot(knoty);
        [n1, n2, n3] = size(w);
        wn     = zeros(n1, n2+1, n3);
        Pxn    = zeros(n1, n2+1, n3);
        Pyn    = zeros(n1, n2+1, n3);
        Pzn    = zeros(n1, n2+1, n3);
        for i = 1:(k - p)
            wn(:,i,:)    = w(:,i,:);
            Pxn(:,i,:)   = Pwx(:,i,:)./wn(:,i,:); 
            Pyn(:,i,:)   = Pwy(:,i,:)./wn(:,i,:); 
            Pzn(:,i,:)   = Pwz(:,i,:)./wn(:,i,:); 
        end
        for i = (k - p + 1):k
            alpha      = (t - knoty(i))/(knoty(i + p) - knoty(i));
            wn(:,i,:)  = (1-alpha)*w(:,i-1,:)+alpha*w(:,i,:);
            Pxn(:,i,:) = (alpha*Pwx(:,i,:)+(1-alpha)*Pwx(:,i-1,:))./wn(:,i,:);
            Pyn(:,i,:) = (alpha*Pwy(:,i,:)+(1-alpha)*Pwy(:,i-1,:))./wn(:,i,:);
            Pzn(:,i,:) = (alpha*Pwz(:,i,:)+(1-alpha)*Pwz(:,i-1,:))./wn(:,i,:);
        end
        for i = (k + 1):(m + 1)
            wn(:,i,:) = w(:,i-1,:);
            Pxn(:,i,:) = Pwx(:,i-1,:)./w(:,i-1,:);
            Pyn(:,i,:) = Pwy(:,i-1,:)./w(:,i-1,:);
            Pzn(:,i,:) = Pwz(:,i-1,:)./w(:,i-1,:);
        end

knotyn = sort([knoty,t]);
end

end