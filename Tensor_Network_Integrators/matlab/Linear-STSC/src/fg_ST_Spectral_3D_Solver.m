function [fgsol,Ctt] = fg_ST_Spectral_3D_Solver(testcase,SN,type,tol,a,b)
%this is a spectral solver for the PDE
% du/dt - K(t,x)*\Deltau + b(t,x)\dot \grad_u + c(t,x)u = rhs(t,x)
% time domain: [0,1]
% space domain: [-1,1]

%% load test case
kfun = testcase.kfun;
bfun = testcase.bfun;
cfun = testcase.cfun;
rhsfn = testcase.rhsfn;
gfn = testcase.gfn;

%%
N1 = SN+1;
%%
% space matrix
[tempDMA,temppts] = compute_derivative_mat_fn(SN,type);

XPt = (b-a)./2.*temppts + (a+b)/2;
DMA = (2/(b-a))*tempDMA;

Lapmat = DMA*DMA; %Laplace matrix

%time matrix
At = 2*tempDMA;
Tvec = (temppts+1)*0.5;
I = eye(N1,N1);

%% create grid in tt
Itt = {ones(1,N1),ones(1,N1),ones(1,N1),ones(1,N1)};
Ctt{1} = cell2core(tt_tensor,{Tvec',ones(1,N1),ones(1,N1),ones(1,N1)});
for ic = 2:4
  temp = Itt;
  temp{ic} = XPt';
  Ctt{ic} = cell2core(tt_tensor,temp);
  % Ctt{ic} = full(Ctt{ic},size(Ctt{ic}));
end
%% build du/dt operatorsin TT format
It = 2:N1;
Ix = 2:N1-1;
% du/dt operator
Btt = matrices_to_tt_matrix_fn({At(It,It),I(Ix,Ix),I(Ix,Ix),I(Ix,Ix)});
B = full(Btt);
Bmaptt = matrices_to_tt_matrix_fn({At(It,:),I(Ix,:),I(Ix,:),I(Ix,:)});
Bmap = full(Bmaptt);
%%
A= B;
Amap = Bmap;

%% %%%%%%%%%%%%%%%% Laplace operator %%%%%%%%%%%
if isa(kfun,'function_handle')
  Laptt = matrices_to_tt_matrix_fn({I(It,It),Lapmat(Ix,Ix),I(Ix,Ix),I(Ix,Ix)}) ...
    + matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),Lapmat(Ix,Ix),I(Ix,Ix)})...
    + matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),I(Ix,Ix),Lapmat(Ix,Ix)});

  Lap = full(Laptt);
  Lapmaptt = matrices_to_tt_matrix_fn({I(It,:),Lapmat(Ix,:),I(Ix,:),I(Ix,:)}) ...
    + matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),Lapmat(Ix,:),I(Ix,:)})...
    + matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),I(Ix,:), Lapmat(Ix,:)});
  Lapmap = full(Lapmaptt);

  % tensorize Kfun
  Kfuntt = amen_cross(Ctt, @(x) cross_fun_nD(x,kfun),tol);
  Kfuntt = tt_get_inner(Kfuntt,{2:N1,2:N1-1,2:N1-1,2:N1-1});
  Koptt = make_tt_to_operator(Kfuntt);

  Kop = full(Koptt);

  %
  % Laptt = amen_mm(Koptt,Laptt,tol);
  % Lapmaptt = amen_mm(Koptt,Lapmaptt,tol);
  Lap = Kop*Lap;
  Lapmap = Kop*Lapmap;
  %
  % Att = round(Att - Laptt,tol);
  % Amaptt = round(Amaptt - Lapmaptt,tol);
  A = A - Lap;
  Amap = Amap - Lapmap;
end

%% Convection operator
if isa(bfun{1},'function_handle')
  boptt = cell(1,3);
  for i = 1:3
    btemptt = amen_cross(Ctt, @(x) cross_fun_nD(x,bfun{i}),tol);
    btemptt = tt_get_inner(btemptt,{2:N1,2:N1-1,2:N1-1,2:N1-1});
    boptt{i} = make_tt_to_operator(btemptt);
    boptt{i} = full(boptt{i});
  end
  %
  Convecmat = DMA;
  Convtt1 = matrices_to_tt_matrix_fn({I(It,It),Convecmat(Ix,Ix),I(Ix,Ix),I(Ix,Ix)});
  Conv1 = boptt{1}*full(Convtt1);
  % Convtt1 = amen_mm(boptt{1},Convtt1,tol);


  Convtt2 = matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),Convecmat(Ix,Ix),I(Ix,Ix)});
  Conv2 = boptt{2}*full(Convtt2);
  % Convtt2 = amen_mm(boptt{2},Convtt2,tol);


  Convtt3 = matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix),I(Ix,Ix), Convecmat(Ix,Ix)});
  Conv3 = boptt{3}*full(Convtt3);
  % Convtt3 = amen_mm(boptt{3},Convtt3,tol);

  Conv = Conv1 + Conv2 + Conv3;

  %
  Convtt1 = matrices_to_tt_matrix_fn({I(It,:),Convecmat(Ix,:),I(Ix,:),I(Ix,:)});
  Conv1 = boptt{1}*full(Convtt1);

  Convtt2 = matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),Convecmat(Ix,:),I(Ix,:)});
  Conv2 = boptt{2}*full(Convtt2);

  Convtt3 = matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),I(Ix,:), Convecmat(Ix,:)});
  Conv3 = boptt{3}*full(Convtt3);

  Convmap = Conv1 + Conv2 + Conv3;

  %
  A = A + Conv;
  Amap = Amap + Convmap;

end

%% Reaction term
if isa(cfun,'function_handle')
  Reacttt = matrices_to_tt_matrix_fn({I(It,It),I(Ix,Ix), I(Ix,Ix),I(Ix,Ix)});
  React = full(Reacttt);
  Reactmaptt = matrices_to_tt_matrix_fn({I(It,:),I(Ix,:),I(Ix,:),I(Ix,:)});
  Reactmap = full(Reactmaptt);

  %% cfun
  cfuntt = amen_cross(Ctt, @(x) cross_fun_nD(x,cfun),tol);
  cfuntt = tt_get_inner(cfuntt,{2:N1,2:N1-1,2:N1-1,2:N1-1});
  coptt = make_tt_to_operator(cfuntt);
  

  %
  React = full(coptt)*React;
  Reactmap = full(coptt)*Reactmap;
  
  %
  A = A + React;
  Amap = Amap + Reactmap;
end

%% compute F interior
ftt = amen_cross(Ctt, @(x) cross_fun_nD(x,rhsfn),tol);
%get the inner of ftt
ftt = tt_get_inner(ftt,{2:N1,2:N1-1,2:N1-1,2:N1-1});
%% compute gBC
% compute BC over the whole domain
G = amen_cross(Ctt, @(x) cross_fun_nD(x,gfn),tol);
Gint = tt_set_zero_boundaries(G,{1,[1,N1],[1,N1],[1,N1]});
Gbc = round(G- Gint,tol);

%% compute the rhs
% RHStt = round(ftt - amen_mv(Amaptt,Gbc,tol),tol);
RHS = full(ftt) - Amap*full(Gbc);
%% solve in TT
fgsol = A\RHS;
%%
