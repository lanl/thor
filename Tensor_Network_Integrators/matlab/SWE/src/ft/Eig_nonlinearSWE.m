function Eig = Eig_nonlinearSWE(Q,d,gacc,~)
    % allocate memory
    Eig = zeros(size(Q));
    %
    h  = Q(1,:);
    hU = Q(2,:);
    hV = Q(3,:);
    %
    U = hU./h;
    V = hV./h;
    %
    a = sqrt(gacc*h);
    %
    temp = sqrt(U.^2+V.^2) + a;
    %
    Eig(1,:) = temp(:);
    Eig(2,:) = temp(:);
    Eig(3,:) = temp(:);
    %
end