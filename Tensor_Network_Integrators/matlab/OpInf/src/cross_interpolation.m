function y = cross_interpolation(x,fun,eps_cr,opt1,opt2)
  %
  if(isa(x,'cell'))
    %
    x0   = opt1;
    verb = opt2;
    %
    nswp = 20;
    %
    y = amen_cross(x, @(xx) fun(xx), eps_cr, "verb", verb, "y0", x0, "nswp", nswp, "trunc_method","svd");
    %
    y = round(y,eps_cr);
    %
  else
    %
    verb = opt1;
    vec  = opt2;
    %
    nswp = 1000;
    %
    y = greedy2_cross(x, @(xx) fun(xx), eps_cr, "verb", verb, "nswp", nswp, "vec", vec);
    %
  end
  %
end