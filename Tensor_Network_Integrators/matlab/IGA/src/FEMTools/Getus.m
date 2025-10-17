function us = Getus(u1, I1, u2, I2)
  N = max([I1(:); I2(:)]);
    us = zeros(N, 1);  % Column vector output

    % Assign values from u1 and u2 to correct positions
    us(I1) = u1;
    us(I2) = u2;

end
