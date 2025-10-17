function [wn, Pxn, Pyn, Pzn, knotzn] = KnotInsertMW3_3D(w, Px, Py, Pz, knotz, t)
 n   = length(t);
 wn = w;
 Pxn = Px;
 Pyn = Py;
 Pzn = Pz;
 knotzn = knotz;
 for l = 1:n
     [wn, Pxn, Pyn, Pzn, knotzn] = KnotInsertW3(wn, Pxn, Pyn, Pzn, knotzn, t(l));
 end
 %% Function for multiple time refinement on 3rd-direction
    function [wn, Pxn, Pyn, Pzn, knotzn] = KnotInsertW3(w, Px, Py, Pz, knotz, t)
        Pwx = Px.*w;
        Pwy = Py.*w;
        Pwz = Pz.*w;
        k      = Findk(knotz, t);
        [m, p] = EvaluateKnot(knotz);
        [n1, n2, n3] = size(w);
        wn     = zeros(n1, n2, n3+1);
        Pxn    = zeros(n1, n2, n3+1);
        Pyn    = zeros(n1, n2, n3+1);
        Pzn    = zeros(n1, n2, n3+1);
        for i = 1:(k - p)
            wn(:,:,i)    = w(:,:,i);
            Pxn(:,:,i)   = Pwx(:,:,i)./wn(:,:,i); 
            Pyn(:,:,i)   = Pwy(:,:,i)./wn(:,:,i); 
            Pzn(:,:,i)   = Pwz(:,:,i)./wn(:,:,i); 
        end
        for i = (k - p + 1):k
            alpha      = (t - knotz(i))/(knotz(i + p) - knotz(i));
            wn(:,:,i)  = (1-alpha)*w(:,:,i-1)+alpha*w(:,:,i);
            Pxn(:,:,i) = (alpha*Pwx(:,:,i)+(1-alpha)*Pwx(:,:,i-1))./wn(:,:,i);
            Pyn(:,:,i) = (alpha*Pwy(:,:,i)+(1-alpha)*Pwy(:,:,i-1))./wn(:,:,i);
            Pzn(:,:,i) = (alpha*Pwz(:,:,i)+(1-alpha)*Pwz(:,:,i-1))./wn(:,:,i);
        end
        for i = (k + 1):(m + 1)
            wn(:,:,i) = w(:,:,i-1);
            Pxn(:,:,i) = Pwx(:,:,i-1)./w(:,:,i-1);
            Pyn(:,:,i) = Pwy(:,:,i-1)./w(:,:,i-1);
            Pzn(:,:,i) = Pwz(:,:,i-1)./w(:,:,i-1);
        end

knotzn = sort([knotz,t]);
end

end