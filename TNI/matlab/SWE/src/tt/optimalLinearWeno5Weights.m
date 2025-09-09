function dr=optimalLinearWeno5Weights(r)
    %
    A = wenokCoeff(r,5);
    a = A(3,:);
    %
    C=wenokCoeff(r,3);
    c=@(i,j) C(i+1,j+1);
    %
    d0=a(5)/c(0,2);
    d2=a(1)/c(2,0);
    d1=(a(4)-c(0,1)*d0)/c(1,2);
    %
    dr=[d0 d1 d2]';
