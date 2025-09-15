function u = solExact(x,y,z,t,gam)
    %
    u = zeros([1,size(x)]);
    %
    n1 = 2*pi; n2 = 2*pi; n3=2*pi; w = 6*pi;
    %
    arg = n1*x + n2*y + n3*z - w*t;
    %
    u(1,:) = sin(arg(:));
    %
end