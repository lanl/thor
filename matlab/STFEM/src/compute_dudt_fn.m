function dudtfn = compute_dudt_fn(uexact)

syms t x y z
u = uexact(t,x,y,z);
dudtfn = matlabFunction(diff(u,t),'Var',{t,x,y,z});

end