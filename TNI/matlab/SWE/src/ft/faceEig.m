function data=faceEig(data,d)
    % get size
    Nx_total = data.Nx_total;
    Ny_total = data.Ny_total;
    % setup index arrays 
    i = [4:Nx_total-3+(d==1)];
    j = [4:Ny_total-3+(d==2)];
    %
    im1 = i - (1==d);
    jm1 = j - (2==d);
    %
    data.eig(:) = 0;
    data.eig(:,i,j,:)  = data.Eig(0.5*(data.QL(:,i,j,:)+data.QR(:,im1,jm1,:)),d,data.gacc,data.H);
    %
end