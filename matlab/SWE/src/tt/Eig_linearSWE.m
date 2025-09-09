function Eig = Eig_linearSWE(Q,data,d)
    %
    Eig = sqrt(data.gacc*data.H)*tt_ones(Q{1}.n);
    %
end