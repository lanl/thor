function dF = full_oderhs(t,F,A,Vfun)
% right hand side function for ODE45
[n,~] = size(A);
F = reshape(F,n,n);
dF = fg_oderhs(t,F,A(1,2),A(1,1),Vfun);
% dF = A*F + F*A' + vt(t);
dF = dF(:);

end