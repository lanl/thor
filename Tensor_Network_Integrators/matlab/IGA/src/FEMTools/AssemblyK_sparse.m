function Kuut = AssemblyK_sparse(Kuue, In1, In2)
    m = max(In1(:));
    n = max(In2(:));
    [n1, n2] = size(In1);
    [~, n3] = size(In2);

    % Preallocate arrays for row, col, and val (over-allocate, trim later)
    max_nz = n1 * n2 * n3;
    row_idx = zeros(max_nz,1);
    col_idx = zeros(max_nz,1);
    val = zeros(max_nz,1);
    count = 0;

    for i = 1:n1
        for j = 1:n2
            for k = 1:n3
                idx1 = In1(i,j);
                idx2 = In2(i,k);
                if idx1 ~= -1 && idx2 ~= -1
                    count = count + 1;
                    row_idx(count) = idx1;
                    col_idx(count) = idx2;
                    val(count) = Kuue(j,k,i);
                end
            end
        end
    end

    % Trim arrays in case we over-allocated
    row_idx = row_idx(1:count);
    col_idx = col_idx(1:count);
    val = val(1:count);

    % Assemble sparse matrix
    Kuut = sparse(row_idx, col_idx, val, m, n);
end
