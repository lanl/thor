function uindex1 = DirichletBC3D_Laplace_Surround(ctpxn, ctpyn, ctpzn, IndexE,...
                                                         bc_type, plot_out)
[n1, n2, n3] = size(ctpxn);
%% Boundary conditions
if bc_type == "1"
    [indexubc1, bc1x, bc1y, bc1z] = FindNodeArrayk1_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc2, bc2x, bc2y, bc2z] = FindNodeArrayk1_3D(n1, ctpxn, ctpyn, ctpzn);
    indexubc = cat(1, indexubc1, indexubc2);
    uindex2 = indexubc(:); 
       
elseif bc_type == "2"
    [indexubc1, bc1x, bc1y, bc1z] = FindNodeArrayk2_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc2, bc2x, bc2y, bc2z] = FindNodeArrayk2_3D(n2, ctpxn, ctpyn, ctpzn);
    indexubc = cat(1, indexubc1, indexubc2);
    uindex2 = indexubc(:); 
    
elseif bc_type == "3"
    [indexubc1, bc1x, bc1y, bc1z] = FindNodeArrayk3_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc2, bc2x, bc2y, bc2z] = FindNodeArrayk3_3D(n3, ctpxn, ctpyn, ctpzn);
    indexubc = cat(1, indexubc1, indexubc2);
    uindex2 = indexubc(:);

elseif bc_type == "12"
    [indexubc1a, bc1ax, bc1ay, bc1az] = FindNodeArrayk1_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc1b, bc1bx, bc1by, bc1bz] = FindNodeArrayk1_3D(n1, ctpxn, ctpyn, ctpzn);
    [indexubc2a, bc2ax, bc2ay, bc2az] = FindNodeArrayk2_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc2b, bc2bx, bc2by, bc2bz] = FindNodeArrayk2_3D(n2, ctpxn, ctpyn, ctpzn);
    indexubc = cat(2, indexubc1a, indexubc1b, indexubc2a, indexubc2b);
    uindex2 = indexubc(:); 
    uindex2 = unique(uindex2);

elseif bc_type == "23"
    [indexubc1a, bc1ax, bc1ay, bc1az] = FindNodeArrayk2_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc1b, bc1bx, bc1by, bc1bz] = FindNodeArrayk2_3D(n2, ctpxn, ctpyn, ctpzn);
    [indexubc2a, bc2ax, bc2ay, bc2az] = FindNodeArrayk3_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc2b, bc2bx, bc2by, bc2bz] = FindNodeArrayk3_3D(n3, ctpxn, ctpyn, ctpzn);
    indexubc = cat(2, indexubc1a, indexubc1b, indexubc2a, indexubc2b);
    uindex2 = indexubc(:); 
    uindex2 = unique(uindex2);

elseif bc_type == "13"
    [indexubc1a, bc1ax, bc1ay, bc1az] = FindNodeArrayk1_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc1b, bc1bx, bc1by, bc1bz] = FindNodeArrayk1_3D(n1, ctpxn, ctpyn, ctpzn);
    [indexubc2a, bc2ax, bc2ay, bc2az] = FindNodeArrayk3_3D(1, ctpxn, ctpyn, ctpzn);
    [indexubc2b, bc2bx, bc2by, bc2bz] = FindNodeArrayk3_3D(n3, ctpxn, ctpyn, ctpzn);
    indexubc = cat(2, indexubc1a, indexubc1b, indexubc2a, indexubc2b);
    uindex2 = indexubc(:); 
    uindex2 = unique(uindex2);
end
% Plot BCs
if plot_out == true
    if bc_type == "1" || bc_type == "2" || bc_type == "3"
        plot3(bc1x, bc1y, bc1z,'r*');
        plot3(bc2x, bc2y, bc2z, 'r*');
    elseif bc_type == "12" || bc_type == "23" || bc_type == "13"
        plot3(bc1ax, bc1ay, bc1az,'r*');
        plot3(bc1bx, bc1by, bc1bz,'r*');
        plot3(bc2ax, bc2ay, bc2az, 'r*');
        plot3(bc2bx, bc2by, bc2bz, 'r*');
    end
end

%% Solve vector
uindex1 = 1:max(IndexE(:));
uindex1(ismember(uindex1, uindex2)) = [];