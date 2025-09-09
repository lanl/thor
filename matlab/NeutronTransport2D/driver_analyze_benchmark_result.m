close all; clear; clc;
run setup2DSlab.m;

%% load result files
experiment_name = 'alpha_problem';
experiment_name = 'keffective';
Ns = [128,256,512];
I = 1:3;
for i = I
  R{i} = load(sprintf('./Results/%s_%d.mat',experiment_name, Ns(i)));
end

%%
gridsize = Ns(I);
numIter = cellfun(@(c) sum(c.Itertime>0),R);
time = cellfun(@(c) sum(c.Itertime),R);
timeperiter = time./numIter;
Hcompress = cellfun(@(c) compress_ratio_tt(c.ttH),R);
Hsize = cellfun(@(c) prod(c.ttH.n.*c.ttH.m)*8,R)/(1024)^7;
err = cellfun(@(c) norm(c.ktt-c.true_k),R);
%err = cellfun(@(c) norm(c.alphatt(end)-c.true_alpha),R);
MemPsi = cellfun(@(c) tt_get_mem(c.ttPsi1), R);

T = table(gridsize',numIter',time',timeperiter',MemPsi', Hsize',Hcompress',err', ...
  'VariableNames',{'Grid Size', 'Number of Iterations', 'Elapsed Time (s)', ...
  'Time per Iter (s)','MemPsi','full grid H in Zetabyte','H compression', 'Eigval Error'})

writetable(T,sprintf('./Results/%s_Table.xlsx',experiment_name))


%% routines
function m = tt_get_mem(var)
a = whos(inputname(1));
m = a.bytes/10^6;
end