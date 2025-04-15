function plot_ranks_error(tollevel, filename)
% Each row of the array of cells is an RK method results 
% Each column is results for a given N

L = load(filename);
data = L.(tollevel); 
[num_methods, num_Nvals] = size(data);
tvals  = 0:0.5:20;      % time samples for solution comparison

% colors
ttcolor = [0.6350 0.0780 0.1840]; %redish
fgcolor = [0 0.4470 0.7410]; %blueish
refcolor = [0 0 0]; %black

for iM = 1:num_methods
    for iN = 1:num_Nvals

        tterror = data{iM,iN}.tterror; 
        fgerror = data{iM,iN}.fgerror;
        maxranks = data{iM,iN}.ttranks;
        refrank = data{iM,iN}.refranks;
        plottitle = [tollevel,'_',data{iM,iN}.mname, '_',num2str(data{iM,iN}.N)];
        fprintf('Plotting %s \n', plottitle)

        figure 
        title(plottitle)
        subplot(1,2,1)
        semilogy(tvals,tterror,'color',ttcolor)
        hold on
        semilogy(tvals,fgerror,'color',fgcolor)
        axis square
        legend('TT-ERK','FG-ERK','Location','best')
        xlabel('Time'), ylabel('Relative Error')
        subplot(1,2,2)
        plot(tvals,maxranks,'color',ttcolor)
        hold on
        plot(tvals,refrank,'Color',refcolor)
        axis square
        legend('TT rank','Reference rank','Location','best')
        xlabel('Time'), ylabel('Max Rank')
        % exportgraphics(gcf,['figs/',plottitle,'.png'],'Resolution',300)
    end
    hold off
    % close all
end


end

