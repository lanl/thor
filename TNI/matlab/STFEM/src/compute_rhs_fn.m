function rhs = compute_rhs_fn(uexact,a)

syms x y;
u = uexact(x,y);
a = a(x,y);

rhs = -divergence(a.*gradient(u));
rhs = matlabFunction(rhs,"Vars",{x,y});

end
