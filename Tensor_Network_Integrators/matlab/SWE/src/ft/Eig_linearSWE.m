function Eig = Eig_linearSWE(Q,d,g,H)
    %
    % allocate memory
    Eig = zeros(size(Q));
    %
    Eig(1,:) = sqrt(g*H);
    Eig(2,:) = sqrt(g*H);
    Eig(3,:) = sqrt(g*H);
    %
end