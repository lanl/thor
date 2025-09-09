function tt=expaxbyct(x,y,t,a,b,c,eps_tt)
    %
    x.core{1} = a*x.core{1};
    y.core{2} = b*y.core{2};
    t         = c*t;
    %
    tt = exp(t)*exptt(x,1).*exptt(y,2);    
    %
    tt = round(tt,eps_tt);
    %
end