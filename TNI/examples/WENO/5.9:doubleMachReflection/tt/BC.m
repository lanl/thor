function Q=BC(Q,data,t,d)
    %
    th     = pi/3;
    gam    = data.gam;
    eps_tt = data.tt.eps;
    eps_cr = data.tt.eps_cr;
    %
    n = Q{1}.n;
    %
    x = data.xg;
    %
    rhoL =  8;
    uL   =  8.25*sin(th);
    vL   = -8.25*cos(th);
    pL   =  116.5; 
    %
    rhoR = 1.4;
    uR   = 0;
    vR   = 0;
    pR   = 1;
    %
    if(d==1)
        %
        % Create TTs for BC enforcement 
        %
        Q1  = cell(1,5); % ghost cells with i=1,2,3
        Q2  = cell(1,5); % ghost cells with i=Nx-2,Nx-1,Nx-1
        %
        invRho = tt_inv(Q{1},eps_tt,eps_cr); 
        %
        % supersonic inflow
        %   
        % start with creating a rank-1 TT to apply BCs on the left
        % First 3 ghost cells are filled with ones and every other cell is set to zero
        %
        tt1 = tt_ones(n);
        %
        cr            = tt1{1};
        cr(:,4:end,:) = 0;
        tt1{1}        = cr;
        %
        % Set the ghost cells values for the inlet
        %
        Q1{1} = rhoL*   tt1;
        Q1{2} = rhoL*uL*tt1;
        Q1{3} = rhoL*vL*tt1;
        Q1{4} = tt_zeros(n);
        Q1{5} = (pL/(gam-1) + 0.5*rhoL*(uL^2+vL^2))*tt1;
        %
        % subsonic outflow (set 4 variables from interior and pressure from the right BC)
        %  
        % these 6 TTs will have ghost cels at the outlet filled by extrapolating from the interior
        % but then interior values are set to zero
        %
        invRhoin = set_outlet_bc(invRho); 
        rhoin    = set_outlet_bc(Q{1});   
        rhoUin   = set_outlet_bc(Q{2});   
        rhoVin   = set_outlet_bc(Q{3});   
        rhoWin   = tt_zeros(n);
        rhoEin   = set_outlet_bc(Q{5});  
        % 
        rhoU2in = round(rhoUin.^2 + rhoVin.^2,eps_tt);
        KE      = round(invRhoin.*rhoU2in,eps_tt);
        %
        pin = round(rhoEin - 0.5*KE,eps_tt)*(gam-1);
        %
        pbc = tt_ones(n); 
        cr  = pbc{1};
        cr(:,1:end-3,:)   = 0;
        cr(:,end-2:end,:) = pR;
        %
        pbc{1} = cr;
        %
        rhoEbc = round((2*pbc-pin)/(gam-1) + 0.5*KE,eps_tt);  
        %
        % Set the ghost cells values for the outlet
        %
        Q2{1} = rhoin;
        Q2{2} = rhoUin;
        Q2{3} = rhoVin;
        Q2{4} = rhoWin;
        Q2{5} = rhoEbc;
        %
        % We need to set the ghost cells for the solution TT to zero before
        % applying the BCs in the x-direction
        %
        Q=zero_ghost(Q,d);
        %
        % This is how BCs are implemented in the x direction in the TT format
        % Q1: inlet BC  --> full tensor form has non-zero entries in Q1(1:3,:,:)
        % Q2: Outlet BC --> full tensor form has non-zero entries in Q2(end-2:end,:,:)
        % Q : solution  --> full tensor form has non-zero entries in Q(4:end-3,:,:)
        %
        for eq=1:5
            Q{eq} = round(Q{eq}+Q1{eq}+Q2{eq},eps_tt);
        end
        %
    elseif(d==2)
        %
        % Create TTs for BC enforcement 
        %
        Q1  = cell(1,5); % for x<1/6 and y=0
        Q2  = cell(1,5); % for x>1/6 and y=0
        Q3  = cell(1,5); % for x<x_shock and y=1
        Q4  = cell(1,5); % for x>x_shock and y=1
        %
        Ly = 1; % y-domain size for this problem
        % 
        % find the x-grid point index with x=1/6 on the bottom boundary
        %
        [~,idx16]=min(abs(x-1/6)); 
        if(x(idx16)>1/6) 
            idx16 = idx16-1;
        end
        %
        % find the shock location on the top boundary and the x-grid point index
        %
        xs  = 10*t/sin(th) + 1/6 + Ly*cot(th);
        [~,idxs]=min(abs(x-xs));  
        if(x(idxs)>xs) 
            idxs = idxs-1;
        end
        %
        % BCs for y=0 & x<1/6 
        %  
        tt1 = tt_zeros(n);
        %
        crx = tt1{1};
        cry = tt1{2};
        crz = tt1{3};
        %
        crx(:,1:idx16,:) = 1; % to set BCs only for x<1/6
        cry(:,1:3,:)     = 1; % to set BCs only for ghost cells with y<0
        crz(:,:,:)       = 1; % to make sure the solution is the uniform in the z-direction
        %
        tt1{1} = crx;
        tt1{2} = cry;
        tt1{3} = crz;
        %
        % Set the ghost cell values for y=0 & x<1/6 
        %
        Q1{1} = rhoL*tt1;
        Q1{2} = rhoL*uL*tt1;
        Q1{3} = rhoL*vL*tt1;
        Q1{4} = tt_zeros(n);
        Q1{5} = (pL/(gam-1) + 0.5*rhoL*(uL^2+vL^2))*tt1;
        %
        % y=0 & x>1/6 (reflective BC)
        %
        % Set the ghost cell values for y=0 & x>1/6 
        %
        Q2{1} =  symmetric_BC(Q{1},idx16); 
        Q2{2} =  symmetric_BC(Q{2},idx16); 
        Q2{3} = -symmetric_BC(Q{3},idx16);  % needed to make sure the BC is reflective for the y-velocity
        Q2{4} =  tt_zeros(n);
        Q2{5} =  symmetric_BC(Q{5},idx16); 
        %  
        % y=Ly & x<x_shock
        %  
        tt1 = tt_zeros(n);
        %
        crx = tt1{1};
        cry = tt1{2};
        crz = tt1{3};
        %
        crx(:,1:idxs,:)    = 1;  % to set BCs only for x<x_shock
        cry(:,end-2:end,:) = 1;  % to set BCs only for ghost cells with y>1
        crz(:,:,:)         = 1;  % to make sure the solution is the uniform in the z-direction
        %
        tt1{1} = crx;
        tt1{2} = cry;
        tt1{3} = crz;
        %
        % Set the ghost cell values for y=1 & x<x_shock
        %
        Q3{1} = rhoL*tt1;
        Q3{2} = rhoL*uL*tt1;
        Q3{3} = rhoL*vL*tt1;
        Q3{4} = tt_zeros(n);
        Q3{5} = (pL/(gam-1) + 0.5*rhoL*(uL^2+vL^2))*tt1;
        %  
        % y=Ly & x>x_shock
        %  
        tt1 = tt_zeros(n);
        %
        crx = tt1{1};
        cry = tt1{2};
        crz = tt1{3};
        %
        crx(:,idxs+1:end,:) = 1;  % to set BCs only for x>x_shock
        cry(:,end-2:end,:)  = 1;  % to set BCs only for ghost cells with y>1
        crz(:,:,:)          = 1;  % to make sure the solution is the uniform in the z-direction
        %
        tt1{1} = crx;
        tt1{2} = cry;
        tt1{3} = crz;
        %
        % Set the ghost cell values for y=1 & x>x_shock
        %
        Q4{1} = rhoR*tt1;
        Q4{2} = rhoR*uR*tt1;
        Q4{3} = rhoR*vR*tt1;
        Q4{4} = tt_zeros(n);
        Q4{5} = (pR/(gam-1) + 0.5*rhoR*(uR^2+vR^2))*tt1;
        %
        % We need to set the ghost cells for the solution TT to zero before
        % applying the BCs in the y-direction
        %
        Q=zero_ghost(Q,d);
        %
        % This is how BCs are implemented in the y direction in the TT format
        % Q1: y=0 & x<1/6     --> full tensor form has non-zero entries in Q1(1:idx16,1:3,:)
        % Q2: y=0 & x>1/6     --> full tensor form has non-zero entries in Q2(idx16+1:end,1:3,:)
        % Q3: y=1 & x<x_shock --> full tensor form has non-zero entries in Q3(1:idxs,end-2:end,:)
        % Q4: y=1 & x>x_shock --> full tensor form has non-zero entries in Q4(1+idxs:end,end-2:end,:,:)
        % Q : solution        --> full tensor form has non-zero entries in Q(:,4:end-3,:)
        %
        for eq=1:5
            Q{eq} = round(Q{eq}+Q1{eq}+Q2{eq}+Q3{eq}+Q4{eq},eps_tt);
        end
        %
    else
        %
        % assuming periodic in z is OK
        %
        Q=BCperiodic(Q,data,t,d);
        %
    end 
    %
end
%
function qbc = set_outlet_bc(qin)
    %
    % Init BC from interior
    %
    qbc = qin;
    %
    % get the TT core for the x-direction
    %
    cr = qbc{1}; 
    %
    % fill ghost cells by extrapolating from the interior
    %
    cr(:,qbc.n(1),:)   = cr(:,qbc.n(1)-3,:);
    cr(:,qbc.n(1)-1,:) = cr(:,qbc.n(1)-3,:);
    cr(:,qbc.n(1)-2,:) = cr(:,qbc.n(1)-3,:);
    %
    % set the interior to zero (needed for the BC implementation in the TT-format)
    %
    cr(:,1:qbc.n(1)-3,:) = 0;
    %
    % This BC only has info about the outlet BC
    %
    qbc{1} = cr;
    %
end
%
function qbc = symmetric_BC(qin,i1)
    %
    % Init BC from interior
    %
    qbc = qin;
    %
    % Set everything to zero for x<1/6 for the x-core
    % 
    cr           = qbc{1};
    cr(:,1:i1,:) = 0;
    qbc{1}       = cr;
    %
    % get the TT core for the y-direction
    %
    cr = qbc{2};
    %
    % Fill ghost cells by applying symmetry
    %
    cr(:,1,:) = cr(:,6,:);
    cr(:,2,:) = cr(:,5,:);
    cr(:,3,:) = cr(:,4,:);
    %
    % set the interior to zero (needed for the BC implementation in the TT-format)
    %
    cr(:,4:end,:) = 0;
    %
    % This BC only has info for x>1/6 and y=0
    %
    qbc{2} = cr;
    %
end
%
function qbc = top_BC(qin)
    %
    % Init BC from interior
    %
    qbc = qin;
    %
    % get the TT core for the y-direction
    %
    cr = qbc{2};
    %
    % fill ghost cells by simple extrapolating from the interior
    %
    cr(:,end  ,:) = cr(:,end-3,:);
    cr(:,end-1,:) = cr(:,end-3,:);
    cr(:,end-2,:) = cr(:,end-3,:);
    %
    % set the interior to zero (needed for the BC implementation in the TT-format)
    %
    cr(:,1:end-3,:) = 0;
    %
    % This BC only has info for y=1
    %
    qbc{2} = cr;
    %
end