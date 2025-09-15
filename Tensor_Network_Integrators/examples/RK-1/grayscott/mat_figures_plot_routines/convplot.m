close all; clear; clc;
%THE MAT FILES THIS DEPENDS ON NEED FIXING
%%Convergence plot
%% Load in data
temp = load('tt_imex_grayscott.mat');
dtvals_imex = temp.dtvals;
errorttend_imex = cellfun(@(c) c.error,temp.R);
lits_imex = cellfun(@(c) c.lits,temp.R);
time_imex = cellfun(@(c) c.time,temp.R);
mname_imex = temp.mname1;
order_imex = 2;

temp = load('tt_dirk_grayscott.mat');
dtvals_dirk = temp.dtvals;
errorttend_dirk = cellfun(@(c) c.error,temp.R);
lits_dirk = cellfun(@(c) c.lits,temp.R);
time_dirk = cellfun(@(c) c.time,temp.R);
mname_dirk = temp.mname;
order_dirk = 2;

%% Set colors and markers 
% colors
imexcolor = [0.6350 0.0780 0.1840]; %redish
dirkcolor = [0 0.4470 0.7410]; %blueish
fgcolor = [0 0 0]; %black

% markers
imexmark = 'square';
dirkmark = 'o';
fgmark = '*';


%% Convergence
p_imex = polyfit(log(dtvals_imex),log(errorttend_imex),1); % get line of best fit
p_dirk = polyfit(log(dtvals_dirk), log(errorttend_dirk),1);

figure('Units', 'inches', 'Position', [1, 1, 16, 9])
subplot(1,2,1)
loglog(dtvals_imex,errorttend_imex,'-','color',imexcolor,'Marker',imexmark,...
    'DisplayName',['TT-IMEX2 (slope =', num2str(p_imex(1),'%.2f'),')'])
hold on
loglog(dtvals_dirk,errorttend_dirk,'-','color',dirkcolor,'Marker',dirkmark,...
    'DisplayName',['TT-DIRK2 (slope =', num2str(p_dirk(1),'%.2f'),')'])
hold off
set ( gca, 'xdir', 'reverse' )
legend('Location','northeast')
xlabel('dt'), ylabel('Rel. Error (Final time)')
axis square
% title(['Convergence'])

%% Runtime
subplot(1,2,2)
loglog(dtvals_imex,time_imex,'-','color',imexcolor,'Marker',imexmark,...
    'DisplayName','TT-IMEX2')
hold on
loglog(dtvals_dirk,time_dirk,'-','color',dirkcolor,'Marker',dirkmark,...
    'DisplayName','TT-DIRK2')
hold off
% legend('Location','best')
xlabel('dt'), ylabel('Runtime')
axis square
% 
% %% Linear iterations
% subplot(1,3,3)
% loglog(dtvals_imex,lits_imex,'-','color',imexcolor,'Marker',imexmark,...
%     'DisplayName','TT-IMEX2')
% hold on
% loglog(dtvals_dirk,lits_dirk,'-','color',dirkcolor,'Marker',dirkmark,...
%     'DisplayName','TT-DIRK2')
% hold off
% % legend
% xlabel('dt'), ylabel('Linear iterations')
% axis square

% Remove extra margins and save as EPS
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
pos = get(gcf, 'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)], 'PaperPosition', [0 0 pos(3) pos(4)]);
% print(gcf, 'figs/appraisalfigs/gs_conv_eff', '-depsc','-r300');
% saveas(gcf,'figs/appraisalfigs/gs_conv_eff','epsc')