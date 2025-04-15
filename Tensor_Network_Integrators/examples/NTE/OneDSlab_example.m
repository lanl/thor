addpath(genpath('../../matlab/NTE/src/'))
addpath(genpath('../../matlab/utils/chebfun/'))
addpath(genpath('../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../matlab/utils/ttfunc/'))

close all; clear; clc;

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
  
  %% TT fixed-point
  ts = datetime;
  [ttPsi1, ttk, ttconvrate] = qtt_1DSlab_EVP_1group_fn(param, tol);
  tt_time = seconds(datetime-ts);
  %%
  fprintf('L = %d \n',L);
  fprintf('ttk = %.5fs \n', ttk)
  fprintf('tt time = %.5fs \n', tt_time)
  fprintf('Psi1 compress ratio = %.5e \n ', compress_ratio_tt(ttPsi1));
end





