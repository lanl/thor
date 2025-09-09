function Q=BC(Q,data,t,d)
    %
    if(d==1) % only apply BCs for once
        % no flow in y
        d = 2;
        n = Q{1}.n(d);
        %
        core1 = Q{1}{d};
        core2 = Q{2}{d};
        core3 = Q{3}{d};
        %
        core1(:,1:3,:) =  core1(:,[6 5 4],:);
        core2(:,1:3,:) =  core2(:,[6 5 4],:);
        core3(:,1:3,:) = -core3(:,[6 5 4],:);
        %
        core1(:,n-2:n,:) =  core1(:,[n-3 n-4 n-5],:);
        core2(:,n-2:n,:) =  core2(:,[n-3 n-4 n-5],:);
        core3(:,n-2:n,:) = -core3(:,[n-3 n-4 n-5],:);
        %
        Q{1}{d} = core1;
        Q{2}{d} = core2;
        Q{3}{d} = core3;
        %
        d=1;
        % periodic in x
        Q=BCperiodic(Q,data,t,d);
        %
    end
    %
end