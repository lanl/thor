function dr=optimalLinearWeno3Weights(r)
    %
    A = wenokCoeff(r,3);
    a = A(2,:);
    %
    C=wenokCoeff(r,2);
    c=@(i,j) C(i+1,j+1);
    %
    d0=a(3)/c(0,1);
    d1=a(1)/c(1,0);
    %
    dr=[d0 d1]';
