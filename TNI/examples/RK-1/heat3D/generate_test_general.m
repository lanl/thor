function testcase = generate_test_general(uexact,a)
testname='heat3D';

%%Include boundary conditions 
syms t x y z 
u = uexact(t,x,y,z);
dudt = diff(u,t); % symbolic derivative wrt t
Lu = a*laplacian(u,[x,y,z]); % symbolic laplacian
f = matlabFunction(simplify(dudt - Lu)); % compute forcing function
%%
testcase.testname = testname;
testcase.exactfn = uexact;
testcase.forcingfn = f;
testcase.boundaryfn = uexact;
testcase.laplacian = matlabFunction(Lu);
% disp(testcase);
end