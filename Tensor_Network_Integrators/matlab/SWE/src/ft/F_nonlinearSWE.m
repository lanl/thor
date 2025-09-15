function F = F_nonlinearSWE(Q,dir,g,~)
    % allocate memory
    F = zeros(size(Q));
    %
    h  = Q(1,:);
    hU = Q(2,:);
    hV = Q(3,:);
    %
    U = hU./h;
    V = hV./h;
    %
    p = 0.5*g*(h.^2);
    %
    if(dir==1)
        %
        F(1,:) = hU(:);
        F(2,:) = hU(:).*U(:) + p(:);
        F(3,:) = hV(:).*U(:);
        %
    elseif(dir==2)
        %
        F(1,:) = hV(:);
        F(2,:) = hU(:).*V(:);
        F(3,:) = hV(:).*V(:) + p(:);
        %
    end
    %
end
