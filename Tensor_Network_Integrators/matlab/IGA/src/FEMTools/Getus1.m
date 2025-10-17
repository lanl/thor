function us = Getus1(u1, I1, IndexE)    
    n1 = length(I1);
    us = zeros(1, max(IndexE(:))); % Initialize the output vector with zeros
    
    % Assign u1 to us based on indices I1
    for i = 1:n1
        us(I1(i)) = u1(i);
    end
end
