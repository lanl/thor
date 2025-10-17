function N3 = NtoN3(N)
    n = length(N);
    N3 = zeros(3, 3 * n);
    N3(1,1:3:end) = N(:);
    N3(2,2:3:end) = N(:);
    N3(3,3:3:end) = N(:);
end