function Index3D = GetIndex3D(Index1D)
[n1, n2] = size(Index1D);
Index3D = zeros(n1,3*n2);
    for i=1:n1
        temp = [3*Index1D(i,:) - 2; 3*Index1D(i,:) - 1; 3*Index1D(i,:)];
        Index3D(i,:) = temp(:);
    end
end