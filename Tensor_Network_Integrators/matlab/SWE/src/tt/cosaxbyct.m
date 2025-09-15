function tt=cosaxbyct(x,y,t,a,b,c,eps_tt)
    %
    x.core{1} = a*x.core{1};
    y.core{2} = b*y.core{2};
    t         = c*t;
    %
    tt = cos(t)*costt(x,1).*costt(y,2) - cos(t)*sintt(x,1).*sintt(y,2) - costt(y,2).*sintt(x,1)*sin(t) - costt(x,1).*sintt(y,2)*sin(t);    
    %
    tt = round(tt,eps_tt);
    %
end