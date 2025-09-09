function rhs = compute_3D_rhs_fn(uexact,a)

syms x y z;
u = uexact(x,y,z);
a = a(x,y,z);

rhs = -divergence(a.*gradient(u));
rhs = matlabFunction(rhs,"Vars",{x,y,z});

end
