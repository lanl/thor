function [fgsol,Ctt] = fg_ST_FD_FD_3D_Solver(testcase,TN,SN,tol,a,b)
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


% [Ax,XPt] = compute_derivative_mat_fn(SN,type);
% [Ax,XPt] = compute_FD_mat_fn(N1,-1,1);
% Lapmat = Ax*Ax; %Laplace matrix

[~,A_G] = Space_FDFD_MAT(N1-1);
dx=(b-a)/(N1-1); % mesh size in x direction
XPt=a:dx:b;
Lapmat = 1/(dx^2)*(-full(A_G));

Ax = compute_FD_mat_fn(N1,a,b);

[At,Tvec] = compute_FD_mat_fn(TN,0,1);

% Identities
I = eye(N1,N1);
Itime = eye(TN,TN);

%% create grid in tt
Itt = {ones(1,TN),ones(1,N1),ones(1,N1),ones(1,N1)};
Ctt{1} = cell2core(tt_tensor,{Tvec,ones(1,N1),ones(1,N1),ones(1,N1)});
for ic = 2:4
  temp = Itt;
  temp{ic} = reshape(XPt,1,N1);
  Ctt{ic} = cell2core(tt_tensor,temp);
end
%% build du/dt operatorsin TT format
It = 2:TN;
Ix = 2:N1-1;
% du/dt operator
Btt = matrices_to_tt_matrix_fn({At(It,It),I(Ix,Ix),I(Ix,Ix),I(Ix,Ix)});
Bmaptt = matrices_to_tt_matrix_fn({At(It,:),I(Ix,:),I(Ix,:),I(Ix,:)});
%%
Att = full(Btt);
Amaptt = full(Bmaptt);

%% %%%%%%%%%%%%%%%% Laplace operator %%%%%%%%%%%
if isa(kfun,'function_handle')
  Laptt = matrices_to_tt_matrix_fn({Itime(It,It),Lapmat(Ix,Ix),I(Ix,Ix),I(Ix,Ix)}) ...
    + matrices_to_tt_matrix_fn({Itime(It,It),I(Ix,Ix),Lapmat(Ix,Ix),I(Ix,Ix)})...
    + matrices_to_tt_matrix_fn({Itime(It,It),I(Ix,Ix),I(Ix,Ix),Lapmat(Ix,Ix)});


  Lapmaptt = matrices_to_tt_matrix_fn({Itime(It,:),Lapmat(Ix,:),I(Ix,:),I(Ix,:)}) ...
    + matrices_to_tt_matrix_fn({Itime(It,:),I(Ix,:),Lapmat(Ix,:),I(Ix,:)})...
    + matrices_to_tt_matrix_fn({Itime(It,:),I(Ix,:),I(Ix,:), Lapmat(Ix,:)});

  % tensorize Kfun
  Kfuntt = amen_cross(Ctt, @(x) cross_fun_nD(x,kfun),tol);
  Kfuntt = tt_get_inner(Kfuntt,{2:TN,2:N1-1,2:N1-1,2:N1-1});
  Koptt = make_tt_to_operator(Kfuntt);

  %
  % Laptt = amen_mm(Koptt,Laptt,tol);
  Laptt = full(Koptt)*full(Laptt);
  % Lapmaptt = amen_mm(Koptt,Lapmaptt,tol);
  Lapmaptt = full(Koptt)*full(Lapmaptt);
  %
  % Att = round(Att - Laptt,tol);
  Att = Att - Laptt;
  % Amaptt = round(Amaptt - Lapmaptt,tol);
  Amaptt = Amaptt - Lapmaptt;
end

%% Convection operator
if isa(bfun{1},'function_handle')
  boptt = cell(1,3);
  for i = 1:3
    btemptt = amen_cross(Ctt, @(x) cross_fun_nD(x,bfun{i}),tol);
    btemptt = tt_get_inner(btemptt,{2:N1,2:N1-1,2:N1-1,2:N1-1});
    boptt{i} = make_tt_to_operator(btemptt);
  end
  %
  Convecmat = Ax;
  Convtt1 = matrices_to_tt_matrix_fn({Itime(It,It),Convecmat(Ix,Ix),I(Ix,Ix),I(Ix,Ix)});
  % Convtt1 = amen_mm(boptt{1},Convtt1,tol);
  Convtt1 = full(boptt{1})*full(Convtt1);

  Convtt2 = matrices_to_tt_matrix_fn({Itime(It,It),I(Ix,Ix),Convecmat(Ix,Ix),I(Ix,Ix)});
  % Convtt2 = amen_mm(boptt{2},Convtt2,tol);
  Convtt2 = full(boptt{2})*full(Convtt2);

  Convtt3 = matrices_to_tt_matrix_fn({Itime(It,It),I(Ix,Ix),I(Ix,Ix), Convecmat(Ix,Ix)});
  % Convtt3 = amen_mm(boptt{3},Convtt3,tol);
  Convtt3 = full(boptt{3})*full(Convtt3);

  Convtt = Convtt1 + Convtt2 + Convtt3;

  %
  Convtt1 = matrices_to_tt_matrix_fn({Itime(It,:),Convecmat(Ix,:),I(Ix,:),I(Ix,:)});
  % Convtt1 = amen_mm(boptt{1},Convtt1,tol);
  Convtt1 = full(boptt{1})*full(Convtt1);

  Convtt2 = matrices_to_tt_matrix_fn({Itime(It,:),I(Ix,:),Convecmat(Ix,:),I(Ix,:)});
  % Convtt2 = amen_mm(boptt{2},Convtt2,tol);
  Convtt2 = full(boptt{2})*full(Convtt2);

  Convtt3 = matrices_to_tt_matrix_fn({Itime(It,:),I(Ix,:),I(Ix,:), Convecmat(Ix,:)});
  % Convtt3 = amen_mm(boptt{3},Convtt3,tol);
  Convtt3 = full(boptt{3})*full(Convtt3);

  Convmaptt = Convtt1 + Convtt2 + Convtt3;

  %
  Att = Att + Convtt;
  Amaptt = Amaptt + Convmaptt;

end

%% Reaction term
if isa(cfun,'function_handle')
  Reacttt = matrices_to_tt_matrix_fn({Itime(It,It),I(Ix,Ix), I(Ix,Ix),I(Ix,Ix)});
  Reactmaptt = matrices_to_tt_matrix_fn({Itime(It,:),I(Ix,:),I(Ix,:),I(Ix,:)});

  %% cfun
  cfuntt = amen_cross(Ctt, @(x) cross_fun_nD(x,cfun),tol);
  cfuntt = tt_get_inner(cfuntt,{2:TN,2:N1-1,2:N1-1,2:N1-1});
  coptt = make_tt_to_operator(cfuntt);

  %
  % Reacttt = amen_mm(coptt,Reacttt,tol);
  Reacttt = full(coptt)*full(Reacttt);
  % Reactmaptt = amen_mm(coptt,Reactmaptt,tol);
  Reactmaptt = full(coptt)*full(Reactmaptt);

  % Att = round(Att+Reacttt,tol);
  Att = Att+Reacttt;
  % Amaptt = round(Amaptt + Reactmaptt,tol);
  Amaptt = Amaptt + Reactmaptt;
end

%% compute F interior
ftt = amen_cross(Ctt, @(x) cross_fun_nD(x,rhsfn),tol);
%get the inner of ftt
ftt = tt_get_inner(ftt,{2:TN,2:N1-1,2:N1-1,2:N1-1});

%% compute gBC
% compute BC over the whole domain
G = amen_cross(Ctt, @(x) cross_fun_nD(x,gfn),tol);
Gint = tt_set_zero_boundaries(G,{1,[1,N1],[1,N1],[1,N1]});
Gbc = round(G- Gint,tol);

%% compute the rhs
% RHStt = round(ftt - amen_mv(Amaptt,Gbc,tol),tol);
RHStt = full(ftt) - Amaptt*full(Gbc);
%% solve in TT
% opts = {'nswp',500};
% [ttsol,amensolvedata] = amen_solve2_for_ST(Att,RHStt,tol,opts);
fgsol = Att\RHStt;
% amensolvedata = 0;
%%
