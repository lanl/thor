function [ctpxn_new, ctpyn_new, ctpzn_new, IndexE_new] = ...
                            GDomain_Index1D(ctpxn, ctpyn, ctpzn, IndexE,...
                                            g_type, plot_out)
[n1, n2, n3] = size(ctpxn);
%% Boundary conditions
if g_type == "1"
    k1 = 1;
    k2 = n1;
    [indexg1, bc1x, bc1y, bc1z] = FindNodeArrayk1_3D(k1, ctpxn, ctpyn, ctpzn);
    [indexg2, bc2x, bc2y, bc2z] = FindNodeArrayk1_3D(k2, ctpxn, ctpyn, ctpzn);
    indexg1 = reshape(indexg1, [1, n2, n3]);
    indexg2 = reshape(indexg2, [1, n2, n3]);
    %Cut the last row
    ctpxn_new = ctpxn(k1:k2-1,:,:);
    ctpyn_new = ctpyn(k1:k2-1,:,:);
    ctpzn_new = ctpzn(k1:k2-1,:,:);
elseif g_type == "2"
    k1 = 1;
    k2 = n2;
    [indexg1, bc1x, bc1y, bc1z] = FindNodeArrayk2_3D(k1, ctpxn, ctpyn, ctpzn);
    [indexg2, bc2x, bc2y, bc2z] = FindNodeArrayk2_3D(k2, ctpxn, ctpyn, ctpzn);
    indexg1 = reshape(indexg1, [n1, 1, n3]);
    indexg2 = reshape(indexg2, [n1, 1, n3]);
    %Cut the last row
    ctpxn_new = ctpxn(:,k1:k2-1,:);
    ctpyn_new = ctpyn(:,k1:k2-1,:);
    ctpzn_new = ctpzn(:,k1:k2-1,:);
elseif g_type == "3"
    [indexg1, bc1x, bc1y, bc1z] = FindNodeArrayk3_3D(k1, ctpxn, ctpyn, ctpzn);
    [indexg2, bc2x, bc2y, bc2z] = FindNodeArrayk3_3D(k2, ctpxn, ctpyn, ctpzn);
    indexg1 = reshape(indexg1, [n1, n2, 1]);
    indexg2 = reshape(indexg2, [n1, n2, 1]);
     %Cut the last row
    ctpxn_new = ctpxn(:,:,k1:k2-1);
    ctpyn_new = ctpyn(:,:,k1:k2-1);
    ctpzn_new = ctpzn(:,:,k1:k2-1);
end
% Plot BCs
if plot_out == true
    plot3(bc1x, bc1y, bc1z,'r*');
    plot3(bc2x, bc2y, bc2z, 'r*');
end

IndexE_new = constrainAndRenumber(IndexE, indexg1, indexg2);


