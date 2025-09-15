%% Info
% Example 5.9 : Double Mach reflection problem
% Developer   : M. Engin Danis
% Notes       : this script can be used to plot 2d numerical solution 

clear all;close all;clc
%% Input
solverType = "tt"; % can be tt or ft 
lvl        = 1;    % need to determine the grid level for the ft result
%% error check
if(solverType~="tt" && solverType~="ft")
    error("solverType can be 'tt' or 'ft'");
end
%
fs  = 19;
%
fig_path = pwd + "/figures";
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end
% load results from full tensor folder
res = load(solverType+"/benchmark-results/Q-lvl="+lvl+".txt");
%
Nx = sqrt(4*size(res,1));
Ny = Nx/4;
%
X   = reshape(res(:,1),Nx,Ny);
Y   = reshape(res(:,2),Nx,Ny);
rho = reshape(res(:,3),Nx,Ny);
%
im = figure;
%
hold on;
%
cmin = 1.5;
cmax = 22.9705;
dc   = (cmax-cmin)/29; % 30 contour levels
%
contourf(X,Y,rho,cmin:dc:cmax);
%
clim([cmin cmax]);
%
colormap(jet(-1));
%
axis equal
%
xlabel("$x$","interpreter","latex");
%
ylabel("$u$","interpreter","latex");
%
title("$\rho$ contours for " + upper(solverType) +" solver","interpreter","latex");
%
set(gca,"FontSize",fs);
%
saveas(im,fig_path+"/solverType="+solverType+"-lvl="+lvl+".png");
%
hold off;
