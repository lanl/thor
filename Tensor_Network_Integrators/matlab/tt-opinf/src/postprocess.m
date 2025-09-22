function res = postprocess(opinflib, filename, opinf_types, isA, isF, tp, th, tint_order, Cdt, skip, eps_tt, ...
                           gamma, g_rom, u, retained_energy, fun_grid1d, func_grid1d_idx)
%
[x,xtt,t,dt,K] = prepareTrainingData(filename,tint_order,th,skip);
%
Neq = size(x,4);
%
res = cell(numel(opinf_types),1);
%
[~,idx]=min(abs(t-tp));
%
tp = t(idx);
%
x_sim = squeeze(x(:,:,:,:,idx));
%
[Nx,Ny,Nz,~,~] = size(x);
%
grid1d = fun_grid1d(Nx,Ny,Nz);
%
grid1d_idx = func_grid1d_idx(Nx,Ny,Nz);
%
fid = fopen("stats.txt","w");
%
for t_idx=1:numel(opinf_types)
  %
  opinf_type = opinf_types{t_idx};
  %
  addpath(opinflib.(opinf_type));
  %
  if(opinf_type=="tt")
    %
    res{t_idx} = opinf(x,u,t,K,isA,isF,tint_order,dt,Cdt*dt,tp,gamma,eps_tt);
    %
    res{t_idx}.xp_ft = full(res{t_idx}.xp,res{t_idx}.xp.n');
    %
  elseif(opinf_type=="qtt")
    %
    xtt2 = tt_reshape(xtt,[prod(xtt.n(1:3)) xtt.n(4) xtt.n(5)]);
    XX   = reshape(x,[prod(xtt.n(1:3)) xtt.n(4) xtt.n(5)]);
    % 
    res{t_idx} = opinf(XX,u,t,K,isA,isF,tint_order,dt,Cdt*dt,tp,gamma,eps_tt);
    %
    res{t_idx}.xp_ft = full(res{t_idx}.xp,res{t_idx}.xp.n');
    %
  elseif(opinf_type=="ft")
    %
    res{t_idx} = opinf(x,u,t,K,isA,isF,tint_order,dt,Cdt*dt,tp,gamma);
    %
    res{t_idx}.xp_ft = reshape(res{t_idx}.xp,[size(x,1) size(x,2) size(x,3) Neq]);
    %
  elseif(opinf_type=="ttrom")
    %
    res{t_idx} = opinf(reshape(x,[numel(x)/size(x,5) size(x,5)]),u,t,K,isA,isF,tint_order,dt,Cdt*dt,tp,g_rom,g_rom,eps_tt);
    %
    res{t_idx}.xp_ft = reshape(res{t_idx}.xp,[size(x,1) size(x,2) size(x,3) Neq]);
    %
  elseif(opinf_type=="rom")
    %
    res{t_idx} = opinf(reshape(x,[numel(x)/size(x,5) size(x,5)]),u,t,K,retained_energy,isA,isF,tint_order,dt,Cdt*dt,tp,g_rom,g_rom);
    %
    res{t_idx}.xp_ft = reshape(res{t_idx}.xp,[size(x,1) size(x,2) size(x,3) Neq]);
    %
  else
    error("Invalid OpInf type: %s",opinf_type)
  end
  %
  xp = squeeze(res{t_idx}.xp_ft(:,:,:,:));
  %
  print_stats(1,res{t_idx},opinf_type,xp,x_sim,t_idx,numel(opinf_types));
  print_stats(fid,res{t_idx},opinf_type,xp,x_sim,t_idx,numel(opinf_types));
  %
  rmpath(opinflib.(opinf_type));
  %
end 
%
fclose(fid);
%
% post processing
%
legends = cell(numel(opinf_types)+1,1);
%
legends{1}="Simulation";
for t_idx = 1:numel(opinf_types)
  legends{t_idx+1} = upper(opinf_types{t_idx})+"-OpInf";
end
%
im = figure('Position', [1000, 100, 800, 600]);
hold on;
plot(grid1d,x_sim(grid1d_idx),"k","LineWidth",2);
%
colors = {'r', 'b', 'g', 'm', 'c', 'k', 'y', [0.5, 0.5, 0.5], [1, 0.5, 0], [0, 0.5, 0.5]};  
markers = {'o', 'x', '^', 's', 'd', '*', '+', 'v', '>', '<', 'p', 'h'};   
%
for t_idx = 1:numel(opinf_types)
  %
  color = colors{mod(t_idx-1, numel(colors)) + 1};
  marker = markers{mod(t_idx-1, numel(markers)) + 1};
  %
  xp = squeeze(res{t_idx}.xp_ft(:,:,:,1));
  %
  plot(grid1d, xp(grid1d_idx), marker, 'Color', color, 'MarkerSize', 12, 'LineWidth', 2);
  %
end
%
xlabel("d");
ylabel("u");
legend(legends);
%
set(gca,"FontSize",18)
grid on
%
saveas(im,"compare-opinf.png");
saveas(im,"compare-opinf.pdf");
%