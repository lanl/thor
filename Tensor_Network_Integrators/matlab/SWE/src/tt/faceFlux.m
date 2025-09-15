function data=faceFlux(data,d)
    %
    % compute + and - numerical conserved variables on the quadrature points at cell faces
    %
    data = data.Recon(data,d);
    %
    % calculate face eigenvalues on the quadrature points at cell faces
    %
    data = faceEig(data,d);
    %
    % calculate the Lax-Friedrichs numerical flux on the quadrature points and perform surface integration
    %
    data = localLF(data,d);
    %
end
