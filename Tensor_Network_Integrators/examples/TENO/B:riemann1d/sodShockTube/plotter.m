%% Info
% Example B   : Sod Shock tube
% Developer   : M. Engin Danis
% Notes       : this script can be used to plot 1d numerical solution 

clear all;close all;clc
%% Input
ft_lvl = 7; % need to determine the grid level for the ft result
tt_lvl = 7; % need to determine the grid level for the tt result
%
fs  = 19;
%
legname = {"FT","TT"};
varname = {"rho","u","p"};
%
fig_path = pwd + "/figures";
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end
% load results from full tensor folder
ft = load("ft/benchmark-results/u-lvl="+ft_lvl+".txt");
% load results from tensor train folder
tt = load("tt/benchmark-results/u-lvl="+tt_lvl+".txt");
%
for var_idx = 1:length(varname)
    %
    im = figure;
    %
    hold on;
    %
    plot(ft(:,1),ft(:,1+var_idx),"b.","LineWidth",2,"MarkerSize",8);
    plot(tt(:,1),tt(:,1+var_idx),"ro","LineWidth",2,"MarkerSize",8);
    %
    xlabel("$x$","interpreter","latex");
    %
    ylabel("$u$","interpreter","latex");
    %
    title(varname{var_idx});
    %
    legend(legname,"Location","northwest","interpreter","latex");
    %
    grid on;
    %
    axis([0 1 0 1]);
    %
    set(gca,"FontSize",fs);
    %
    saveas(im,fig_path+"/"+varname{var_idx}+"-lvl="+ft_lvl+".png");
    %
    hold off;
end
