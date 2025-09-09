close all; clear; clc;
run setup1DSlab.m;

%%
Nl = [2,4,8,16,32];
R = cell(1,numel(Nl));
for ix = 1:numel(Nl)
  L = Nl(ix);
  fprintf('********** L = %d ************** \n', L);
  % Slab Dimensions
  x0 = 0.0; x1 =2*1.853722;

  ibl = 0;
  ibr = 0;

  % Cross Section Data
  sigt = 0.32640;
  sigs = 0.225216;
  nusigf = 3.24*0.081600;
  chi = 1.0;
  vel = 1.0;
  G = length( sigt );
  % Problem Parameters
  M = 1023;
  ords = L;

  % Fission Neutron Group PDF

  % Quadrature Information
  [ mu, wgt ] = AngularQuad1DSlab ( L, -1, 1 );

  %Source
  Q  = ones( G*(M+1)*ords, 1 );
  %Include Fission Source
  nosigf = 0;
  nosigs = 0;


  %% setup parameters
  param.a = 0.0;
  param.b = 2*1.853722;
  % param.alpha = 0; % no second term
  % scatterting parameters
  param.sigma_s = 0.225216;
  % fission parameters
  param.nusigma_f = 3.24*0.081600;
  param.chi = 1.0;

  param.sigma =[0.32640];

  param.nmu = L; % index l
  param.nx = M+1; % index i - number of EDGES

  tol = 1e-6; %both convergence and tt
  param.fixed_point_tol = 1e-6;
  param.niter = 200;
  %% full-grid fixed-point
  ts = datetime;
  [Psi1, fgk, convrate] = fg_1DSlab_EVP_fn(param, tol);
  fg2_time = seconds(datetime-ts);
  %% solver
  ts = datetime;
  [ psi_e, k] = Solve1DSlab ( sigt, sigs, nusigf, chi, vel,...
    mu, wgt, Q, x0, x1, M, ibl, ibr, nosigf, nosigs );
  fg_time = seconds(datetime-ts);
  %   fprintf('fixed point in TT - ktt = %.5e \n', k);
  %% TT fixed-point
  ts = datetime;
  [ttPsi1, ttk, ttconvrate] = qtt_1DSlab_EVP_1group_fn(param, tol);
  tt_time = seconds(datetime-ts);

  %%
  fprintf('||k-ktt|| = %.5e \n', abs(k-ttk));
  temp = sort(full(ttPsi1));
  psi = sort(psi_e);
  fgpsi = sort(Psi1(:));
  fprintf('||Psitt-Psi_e||/||psi|| = %.5e \n', norm(psi-temp)/norm(psi));
  fprintf('||fgPsi-Psi_e||/||psi|| = %.5e \n', norm(psi-fgpsi)/norm(psi));
  fprintf('||fgPsi-ttPsi||/||fgPsi|| = %.5e \n', norm(temp-fgpsi)/norm(fgpsi));
  fprintf('fg time = %.5fs \n', fg_time)
  fprintf('tt time = %.5fs \n', tt_time)
  fprintf('Psi1 compress ratio = %.5e \n ', compress_ratio_tt(ttPsi1));
  %% export results
  R{ix}.psi = psi_e;
  R{ix}.time = fg_time;
  R{ix}.fgPsi = Psi1;
  R{ix}.fgk = fgk;
  R{ix}.k = k;
  R{ix}.fg_time = fg2_time;
  R{ix}.ttPsi1 = ttPsi1;
  R{ix}.ttk = ttk;
  R{ix}.tt_time = tt_time;
  R{ix}.tol = tol;
  R{ix}.param = param;

  save('./Results/1Dslab_Pu_239_results_darwin.mat','R','Nl','-v7.3');
  fprintf('The result is saved for L = %d \n',L);
end





