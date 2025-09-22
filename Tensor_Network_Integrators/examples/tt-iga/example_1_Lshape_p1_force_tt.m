close all; clear; clc;

addpath(genpath('../../matlab/utils/tt-toolbox/'));
addpath(genpath('../../matlab/utils/ttfunc/'));
addpath(genpath('../../matlab/tt-iga/src/'));

%% 1) Define list of nm values
nm_list = [2, 4, 8, 16];
nRun    = numel(nm_list);

%% 2) Constant parameters
LL        = 1;
p1 = 1.0; p2 = 1.0; p3 = 1.0;
u_in      = 0.0;
u_out     = 0.0;
bc_type   = "12";
ffunction = "sinpxsinpy";
tt_tol    = 1e-12;

outFile   = 'Lshape_p1_tt.csv';
if isfile(outFile)
  delete(outFile);
end
%% 3) Loop and save per iteration
for idx = 1:nRun
    fprintf("==================================================\n");
    nm = nm_list(idx);
    % 3a) Build geometry / mesh sizes
    [wn, xN, yN, zN, k1, k2, k3] = LShape111v2(LL, nm, nm, nm);
    [n1, n2, n3] = size(wn);
    wn = ones(size(wn));    % B-spline weights
    fprintf("Lshape p1 tt - row %d/%d: \n", idx, nRun);
    fprintf("Size %d x %d x %d: \n", n1, n2, n3);
    % 3b) IGA connectivity
    [xV, yV, zV, IndexE, xE, yE, zE, we, m3D] = ...
        IGAinfo(xN, yN, zN, wn, k1, k2, k3);

    % 3c) Boundary-condition indices
    uidx = DirichletBC3D_Laplace_Surround(xN, yN, zN, IndexE, bc_type, false);

    %% 3d) TT solve
    tstart = tic;
    [u1tt, ustt, Kcm, fcm, ucm] = LinearSolveLaplaceTTf_3D_cross(k1, k2, k3,...
                                                 xN, yN, zN,...
                                                 u_in, u_out, bc_type,...
                                                 ffunction, tt_tol);
    tt_time = toc(tstart);
    
    % Compute TT L2-error
    %%
    uf   = full(ustt);
    ueTT = uf(IndexE);
    % save('ueTT.mat', 'ueTT', '-v7.3')
    tic
    fprintf("Lshape tt : Compute L2 norm \n");
    L2tt = CalH1L2_Lshape(ueTT, xE, yE, zE, ...
                    k1, k2, k3, ...
                    p1, p2, p3, we, m3D);
    %if n1*n2*n3 > 1e6
    %    L2tt = CalH1L2_Lshape_par(ueTT, xE, yE, zE, ...
    %               k1, k2, k3, ...
    %               p1, p2, p3, we, m3D);
    %else
    %    L2tt = CalH1L2_Lshape(ueTT, xE, yE, zE, ...
    %                k1, k2, k3, ...
    %                p1, p2, p3, we, m3D);
    %end
    fprintf("Lshape p1 tt - L2 norm = %.2e\n", L2tt);
    toc

    %%
%    %% 3e) Full-grid solve
%     tstart = tic;
%     [ueL, ~, ~, ~, ~] = LinearSolveLaplacef_3D(...
%                             xE, yE, zE, ...
%                             k1, k2, k3, ...
%                             p1, p2, p3, ...
%                             we, m3D, IndexE, uidx, ffunction);
%     fg_time = toc(tstart);
% 
%     % Compute FG L2-error
%     [L2fg, ~] = CalH1L2_Lshape(...
%                     ueL, xE, yE, zE, ...
%                     k1, k2, k3, ...
%                     p1, p2, p3, we, m3D);

    %% 3f) Pack into 1-row table and write/append
    T = table(n1, n2, n3, Kcm, fcm, ucm, tt_time, L2tt, 'VariableNames',{ ...
             'n1','n2','n3','K_compression','f_compression',...
             'u_compression','TT_solve_time','L2_TT'} );

    if idx == 1
        % Write header + first row
        writetable(T, outFile);
    else
        % Append just the data row
        writetable(T, outFile, ...
            'WriteMode','append', ...
            'WriteVariableNames',false);
    end

    fprintf('Lshape tt completed size (%d x %d x %d) in %d seconds \n', n1, n2, n3, tt_time);
    fprintf("==============================================================\n");
    %% 3g) Release large temporaries before next iteration
    clearvars -except nm_list nRun LL p1 p2 p3 u_in u_out bc_type ffunction tt_tol outFile idx

end

fprintf('Lshape tt - All runs done. Results are in "%s".\n', outFile);
