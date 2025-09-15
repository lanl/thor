function tt=shiftDiff(tt,k,dir)
    % 
    % 1. This functions returns (I-S_dir)tt
    % 2. Apply the shift operator, S_dir*tt, in a given dimension 
    %    and in the direction of dir
    % 3. Compute the difference, tt = tt - S_dir*tt, 
    %
    idx = (3:tt.n(k)-2); 
    %
    corek = tt{k};
    %
    corek(:,idx,:) = corek(:,idx,:) - corek(:,idx+dir,:);
    %
    tt{k} = corek;
    %
end