function Dfn = compute_3D_reaction_fn(uexact,c)

syms t x y z;
u = uexact(t,x,y,z);
c = c(x,y,z);

rhs = c.*u;
Dfn = matlabFunction(rhs,"Vars",{t,x,y,z});

end
