function data=faceFlux(data,d)
    %
    % compute + and - numerical fluxes for finite difference scheme
    %
    % calculate split fluxes using local Lax-Friedrichs scheme
    data = localLF(data,d);
    % reconstruct split fluxes at cell faces
    data = recon(data,d);
    %
end
  