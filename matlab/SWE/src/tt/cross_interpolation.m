function y=cross_interpolation(x,fun,eps_cr,x0,verb)
    %
    nswp = 20;
    %
    y = amen_cross_swe(x,  @(x) fun(x), eps_cr, "verb", verb, "y0", x0,"nswp",nswp,"trunc_method","svd");
    %
    y = round(y,eps_cr);
