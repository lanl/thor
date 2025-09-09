function Q=BC(Q,data,t,d)
    %
    if(d==1)
        %
        % BCs at x=0 and x=0.025
        %
        Q{1} = BC_x(Q{1},d,1);
        Q{2} = BC_x(Q{2},d,-1);
        Q{3} = BC_x(Q{3},d,1);
        Q{4} = BC_x(Q{4},d,1);
        Q{5} = BC_x(Q{5},d,1);
        %
    elseif(d==2)
        %
        gam = data.gam;
        %
        % Set ghost cell values to zero in the y-direction
        %
        Q = zero_ghost(Q,d);
        %
        % crate a TT to enforce BCs on the ghost cells
        %
        Qbc = cell(1,5);
        %
        % create a temporary rank-1 TT with ones in the x & z cores and zeros in the y core
        %
        temp    = tt_ones(Q{1}.n);
        temp{2} = zeros(size(temp{2}));
        %
        % start assembling the pressure and density values of the ghost cells in the y-direction
        %
        pbc   = temp;
        rhobc = temp;
        %
        % get the y-cores of density and pressure, which are filled with zeros for now
        %
        crho = rhobc{2};
        cp   = pbc{2};
        %
        % BC at y=0
        %
        crho(:,1:3,:) = 2;
        cp(:,1:3,:)   = 1;
        %
        % BC at y=1
        %
        crho(:,end-2:end,:) = 1;
        cp(:,end-2:end,:)  = 2.5;
        %
        % set the y-cores of density and pressure with the updated ghost cell values
        %
        rhobc{2} = crho;
        pbc{2}   = cp;
        %
        % update the ghost cells of the BC-TT while keeping the interior cell values zero
        %
        Qbc{1} = rhobc;       % density BC
        Qbc{2} = temp;        % velocity BC (zero) 
        Qbc{3} = temp;        % velocity BC (zero) 
        Qbc{4} = temp;        % velocity BC (zero) 
        Qbc{5} = pbc/(gam-1); % total energy BC 
        %
        % add ghost cell values 
        %
        for eq=1:5
            Q{eq} = Q{eq}+Qbc{eq};
        end
        %
    elseif(d==3)
        %
        % assuming periodic in z is OK
        %
        Q = BCperiodic(Q,data,t,d);
        %
    end
    %
end
%
function Q = BC_x(Q,d,coeff)
    %
    core = Q{d};
    n    = Q.n(d);
    %
    for i=1:3
        %
        core(:,i,:)     = coeff*core(:,6-i+1,:);
        core(:,n-i+1,:) = coeff*core(:,n-6+i,:);
        %
    end
    %
    Q{d} = core;
    %
end