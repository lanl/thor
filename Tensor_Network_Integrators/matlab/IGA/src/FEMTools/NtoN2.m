function N2 = NtoN2(N)
    N = N(:)';
    n = length(N);

    N2 = zeros(2, 2 * n);

    N2(1, 1:2:end) = N;
    N2(2, 2:2:end) = N;
end