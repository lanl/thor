function generate_test_case_fn(testname, eqntype)


syms t x y z
switch testname
  case 'rank1' %rank-1 seprated 
    %true solution
    u(t,x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z)*sin(pi*t);
  case 'rankr'
    u(t,x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z)*sin(pi*t);
    for k = 2:3
      u(t,x,y,z) = u(t,x,y,z) + sin(k*pi*x)*sin(k*pi*y)*sin(k*pi*z)*sin(k*pi*t);
    end
    %%%%%%%%%%%%%%%%%
  case 'sinpoly2'
    % u(t,x,y,z) = sin(2*pi*(t+x+y+z));
    u(t,x,y,z) = sin(pi.*(t+x+y+z).*(2.0./1.1e+1))+ 1/16*(t+x+y+z).^2;
  case 'poly10'
    u(t,x,y,z) = (x+y+z+t).^10;
  case 'sinexp'
    u(t,x,y,z) = sin(10*pi.*(t+x+y+z))+ 1/16*exp(t+x+y+z);
  case 'abs'
    u(t,x,y,z) = sin(pi*x).*sin(pi*y).*sin(pi*z) + x.^2.*abs(x);
  case 'exp'
    u(t,x,y,z) = exp(-pi*(t+x+y+z));
end

%%

switch eqntype
  case 'Laplace'
    dudt = diff(u,t);
    Lu = laplacian(u,[x,y,z]);
    kfun = @(t,x,y,z) 1;
    bfun = {0,0,0};
    cfun = 0;
    testcase.rhsfn = matlabFunction(dudt - Lu);
  case 'consACDR'
    kfun = @(t,x,y,z) 1;
    bfun1 = @(t,x,y,z) 1;
    bfun2 = @(t,x,y,z) 1;
    bfun3 = @(t,x,y,z) 1;
    bfun = {bfun1,bfun2,bfun3};
    cfun = @(t,x,y,z) 1;
    testcase.rhsfn = generate_ACDR_test_case_fn(kfun,bfun,cfun,u);
  case 'fullACDR'
    kfun = @(t,x,y,z) exp(-t.^2);
    bfun1 = @(t,x,y,z) sin(x);
    bfun2 = @(t,x,y,z) exp(-(x+z).^2);
    bfun3 = @(t,x,y,z) cos(y*z);
    bfun = {bfun1,bfun2,bfun3};
    cfun = @(t,x,y,z) cos(2*pi*(x+y+z+t));
    testcase.rhsfn = generate_ACDR_test_case_fn(kfun,bfun,cfun,u);
  case 'testACDR'
    kfun = @(t,x,y,z) exp(-t.^2);
    bfun1 = @(t,x,y,z) sin(2*pi*x);
    bfun2 = @(t,x,y,z) cos(2*pi*y);
    bfun3 = @(t,x,y,z) sin(2*pi*z);
    bfun = {bfun1,bfun2,bfun3};
    cfun = @(t,x,y,z) cos(2*pi*(x+y+z+t));
    testcase.rhsfn = generate_ACDR_test_case_fn(kfun,bfun,cfun,u);
end
%%
testcase.exactfn = matlabFunction(u);
testcase.kfun = kfun;
testcase.bfun = bfun;
testcase.cfun = cfun;
testcase.gfn = testcase.exactfn;
%% save
disp(testcase)
save(sprintf('./testcases/%s-%s.mat',eqntype,testname),'testcase')
