function m2D = macro2DT(knot1, knot2)
    k1 = ReduceKnotFull(knot1);
    k2 = ReduceKnotFull(knot2);
    
    % Get the sizes
    len1 = length(k1) - 1;
    len2 = length(k2) - 1;
    
    % Initialize the 2D cell array
    m2D = zeros(len1, len2, 4);
    
    % Loop through and populate the cell array
    for i = 1:len1
        for j = 1:len2
            m2D(i, j, :) = [k1(i) k1(i+1) k2(j) k2(j+1)];
        end
    end
end
