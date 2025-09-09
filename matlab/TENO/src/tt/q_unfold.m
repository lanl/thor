function Q = q_unfold(Qc,eps)
    %
    d = Qc.d;
    r = Qc.r;
    n = Qc.n;
    %
    Neq = r(d+1);
    %
    if(Neq==1)
        Q{1} = Qc;
        return;
    end
    %
    r(d+1) = 1;
    %
    Q = cell(1,Neq);
    %
    for eq=1:Neq
        %
        Q{eq} = create_tt(n,r,d);
        %
        Q{eq}{1} = Qc{1};
        Q{eq}{2} = Qc{2};
        %
        core = Q{eq}{3};
        %
        core(:,:,:) = Qc{3}(:,:,eq);
        %
        Q{eq}{3} = core;
        %
        Q{eq} = round(Q{eq},eps);
        %
    end
    %
end
