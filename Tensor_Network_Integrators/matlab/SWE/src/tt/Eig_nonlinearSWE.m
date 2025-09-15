function Eig = Eig_nonlinearSWE(Q,data,d)
    %
    if(d<3) % used for numerical flux computations
        %
        % check if initial guess has been set
        %
        [~,idx] = getMaxRank(Q);
        %
        eig0 = Q{idx}; % this is only an initial guess to facilitate the cross interpolation
        %
        Eig = cross_interpolation(Q, @(x) funEig(x,d,data.gacc), data.tt.eps_cr, eig0, 0);
        %
    else % used to determine the time step
        %
        [~,idx] = getMaxRank(data.tt.Q);
        %
        eig0 = data.tt.Q{idx}; % this is only an initial guess to facilitate the cross interpolation
        %
        Eig = cross_interpolation(Q, @(x) funEig(x,d,data.gacc), data.tt.eps_cr, eig0, 0);
        %
    end
    %
end
%
function Eig = funEig(Q,d,g)
    %
    h  = Q(:,1);
    hU = Q(:,2);
    hV = Q(:,3);
    %
    U = hU./h;
    V = hV./h;
    %
    a = sqrt(g*h);
    %
    if(d==1)
        Eig = abs(U) + a;
    elseif(d==2)
        Eig = abs(V) + a;
    else
        Eig = sqrt(U.^2 + V.^2) + a;
    end
    %
end