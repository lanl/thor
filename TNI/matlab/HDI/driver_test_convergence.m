%%
% test convergence

%% %%%%%
close all; clear; clc;
run setup.m;
% run setup_advanpix_drac.m;
run setup_advanpix_temp.m;

%%
mpToggle = 0; %use mp or not 

norder = 2;
Ns = [2];
d = 4;
q = 3;
T = zeros(1,numel(Ns));

if mpToggle
  tol = 1e-30;
  Errors = mp(zeros(1,numel(Ns)));
  for i = 1:numel(Ns)
    fprintf('\n\n --- Computing wit ncell = %d -------- \n\n',Ns(i));
    
    tic
    Errors(i) = tt_GL_mp_fn(d, q, norder, Ns(i), tol);
    T(i) = toc;
    fprintf('\n Execution time for ncell = %d : %.2f seconds \n', Ns(i), T(i))
  end
else % no mp
  tol = 1e-15;
  Errors = zeros(1,numel(Ns));
  for i = 1:numel(Ns)
    fprintf('\n\n --- Computing wit ncell = %d -------- \n\n',Ns(i));
    tic
    Errors(i) = tt_GL_fn_poly(d, q, norder, Ns(i), tol);
    T(i) = toc;
    fprintf('\n Execution time for ncell = %d : %.2f seconds \n', Ns(i), T(i))
  end
end



%%
clc;
fprintf('-------- Convergence rate, d = %d, q = %d, GL order = %d ------------\n',d,q,norder)
fprintf('n = %d, err = %.2e, time = %.2fs, conv = NaN \n',Ns(1), Errors(1), T(1))
for i = 2:numel(Ns)
  conv = log(Errors(i)/Errors(i-1))/log(Ns(i-1)/Ns(i));
  fprintf('n = %d, err = %.2e, time = %.2fs, conv = %.5f \n',Ns(i), Errors(i), T(i), conv)
end