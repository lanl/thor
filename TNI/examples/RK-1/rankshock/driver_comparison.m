Methods  = {'Bogacki-Shampine-ERK','Fehlberg-ERK','Zonneveld-4-3-ERK'};
Nvals = [100:100:500,1000];
loose_tol = [1e-4,1e-9];
medium_tol = [1e-5,1e-10]; 
tight_tol = [1e-6,1e-11];
dtvals = [0.1,1e-6,0.5]; % initial, min, max

% Loose tol
for iM = 1:length(Methods)
    for iN = 1:length(Nvals)
        N = Nvals(iN);
        fprintf('N = %d.  Method = %s. \n',N,Methods{iM});
        loose{iM,iN} = rankshock_comparison(Methods{iM},dtvals,N,loose_tol(1),loose_tol(2));
    end
    save('mat_figures_plot_routines/loose_tol_results.mat','loose');
end

% Medium tol
for iM = 1:length(Methods)
    for iN = 1:length(Nvals)
        N = Nvals(iN);
        fprintf('N = %d.  Method = %s. \n',N,Methods{iM});
        medium{iM,iN} = rankshock_comparison(Methods{iM},dtvals,N,medium_tol(1),medium_tol(2));
    end
    save('mat_figures_plot_routines/medium_tol_results.mat','medium');
end

% Tight tol
for iM = 1:length(Methods)
    for iN = 1:length(Nvals)
        N = Nvals(iN);
        fprintf('N = %d.  Method = %s. \n',N,Methods{iM});
        tight{iM,iN} = rankshock_comparison(Methods{iM},dtvals,N,tight_tol(1),tight_tol(2));
    end
    save('mat_figures_plot_routines/tight_tol_results.mat','tight');
end