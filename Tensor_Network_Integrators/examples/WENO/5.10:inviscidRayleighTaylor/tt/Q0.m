function data=Q0(data)
    %
    gam = data.gam;
    %
    x = data.xg;
    y = data.yg;
    %
    n = data.tt.Q{1}.n;
    %
    % determine the indices below and above the initial discontinuity
    %
    idxL = y<0.5;
    idxR = y>=0.5;
    % 
    % calculate rank-1 cosine term
    %
    cosX    = tt_ones(n);
    cosX{1} = cos(8*pi*x);
    %
    % initialize density, v-velocity and pressure (not ready now)
    %
    rho  = tt_ones(n);
    v    = tt_ones(n);
    p    = tt_ones(n);
    %
    % get the y-cores of density and pressure, which are filled with ones for now
    %
    crho = rho{2};
    cp   = p{2};
    %
    % set values below the discontinuity
    %
    crho(:,idxL,:) = 2;
    cp(:,idxL,:)   = 2*y(idxL)+1;
    %
    % set values above the discontinuity
    %
    crho(:,idxR,:) = 1;
    cp(:,idxR,:)   = y(idxR)+1.5;
    %
    % finalize the initial density and pressure values (since this is a rank-1 IC)
    %
    rho{2} = crho;
    p{2}   = cp;
    %
    % calculate the y-velocity by manipulating the y-core of v-velocity (since this is a rank-1 IC)
    %
    v{2} = -0.025*sqrt(gam*p{2}./rho{2});
    %
    % complete the y-velocity by multiplying with the cosine term
    %
    v = round(v.*cosX,data.tt.eps);
    %
    % set the IC for the conserved variables
    %
    data.tt.Q{1} = rho;
    data.tt.Q{2} = tt_zeros(n);
    data.tt.Q{3} = round(rho.*v,data.tt.eps);
    data.tt.Q{4} = tt_zeros(n);
    data.tt.Q{5} = round(p/(gam-1) + round(0.5*data.tt.Q{3}.*v,data.tt.eps),data.tt.eps);
    %
end
