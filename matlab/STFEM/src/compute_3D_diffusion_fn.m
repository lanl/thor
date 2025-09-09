function Dfn = compute_3D_diffusion_fn(uexact,a)

syms t x y z;
u = uexact(t,x,y,z);
a = a(x,y,z);

gradterm = a.*gradient(u,[x,y,z]);
rhs = -divergence(gradterm,[x,y,z]);
Dfn = matlabFunction(rhs,"Vars",{t,x,y,z});

end
