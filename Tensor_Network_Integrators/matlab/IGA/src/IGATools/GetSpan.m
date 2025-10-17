function idx = GetSpan(knot,zita)

[n, p] =EvaluateKnot(knot);

    if (zita == knot(n + 1))
        idx = n;
        return
    end
    low = p + 1; 
    high = n + 1;
    mid = floor((low + high)/2);
    while(zita < knot(mid) ||...
          zita >= knot(mid + 1))
       if(zita < knot(mid))
            high = mid;
       else
            low = mid;
       end
        mid = floor((low + high) / 2);
    end
    idx = mid;

end

