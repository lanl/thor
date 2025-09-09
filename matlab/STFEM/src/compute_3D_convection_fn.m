function Cvecfn = compute_3D_convection_fn(uexact,bfn)

syms t x y z;
u = uexact(t,x,y,z);
bx = bfn{1}(x,y,z);
by = bfn{2}(x,y,z);
bz = bfn{3}(x,y,z);

gradterm = [bx;by;bz].*gradient(u,[x,y,z]);


Cvecfn = matlabFunction(sum(gradterm),"Vars",{t,x,y,z});

end
