function [uindex1, uindex2, ubc] = DirichletBC3D_Laplace(ctpxn, ctpyn, ctpzn, IndexE,...
                                                         bc_type, val1, val2, plot_out)
[n1, n2, n3] = size(ctpxn);
%% Boundary conditions
if bc_type == "1"
    [indexubc1, bc1x, bc1y, bc1z] = FindNodeArrayk1_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc2, bc2x, bc2y, bc2z] = FindNodeArrayk1_3D(n1, ctpxn, ctpyn, ctpzn);
    indexubc1 = reshape(indexubc1, [1, n2, n3]);
    indexubc2 = reshape(indexubc2, [1, n2, n3]);
    indexubc = cat(1, indexubc1, indexubc2);
    uindex2 = indexubc(:); 
    ubc1 = val1*ones(1, n2, n3);
    ubc2 = val2*ones(1, n2, n3);
    ubcten = cat(1, ubc1, ubc2);
elseif bc_type == "2"
    [indexubc1, bc1x, bc1y, bc1z] = FindNodeArrayk2_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc2, bc2x, bc2y, bc2z] = FindNodeArrayk2_3D(n2, ctpxn, ctpyn, ctpzn);
    indexubc1 = reshape(indexubc1, [n1, 1, n3]);
    indexubc2 = reshape(indexubc2, [n1, 1, n3]);
    indexubc = cat(1, indexubc1, indexubc2);
    uindex2 = indexubc(:); 
    ubc1 = val1*ones(n1, 1, n3);
    ubc2 = val2*ones(n1, 1, n3);
    ubcten = cat(1, ubc1, ubc2);
elseif bc_type == "3"
    [indexubc1, bc1x, bc1y, bc1z] = FindNodeArrayk3_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc2, bc2x, bc2y, bc2z] = FindNodeArrayk3_3D(n3, ctpxn, ctpyn, ctpzn);
    indexubc1 = reshape(indexubc1, [n1, n2, 1]);
    indexubc2 = reshape(indexubc2, [n1, n2, 1]);
    indexubc = cat(2, indexubc1, indexubc2);
    uindex2 = indexubc(:);
    ubc1 = val1*ones(n1, n2, 1);
    ubc2 = val2*ones(n1, n2, 1);
    ubcten = cat(2, ubc1, ubc2);
end
% Plot BCs
if plot_out == true
    plot3(bc1x, bc1y, bc1z,'r*');
    plot3(bc2x, bc2y, bc2z, 'r*');
end

%% BC vector
ubc = ubcten(:)';

%% Solve vector
uindex1 = 1:max(IndexE(:));
uindex1(ismember(uindex1, uindex2)) = [];
