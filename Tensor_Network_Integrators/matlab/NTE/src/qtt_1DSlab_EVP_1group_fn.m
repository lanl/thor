function [qttPsi1, k, convrate, Itertime] = qtt_1DSlab_EVP_1group_fn(param, tol)

%% Load parameters
a = param.a;
b = param.b;

sigma_s = param.sigma_s;
nusigma_f = param.nusigma_f;
chi = param.chi;
sigma = param.sigma;
nmu = param.nmu;
nx = param.nx;
nE = numel(sigma); % index g
dx = (b-a)/(nx-1); %spacial step
fixed_point_tol = param.fixed_point_tol;
%% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]);
mu = flip(mu');
wgt = flip(wgt'/2);

%% forming the system in the tensor format
% positive mu part
Lposmumat = 1/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));
Lposgmat = diag(ones(nx,1)) + diag(ones(nx-1,1),-1);

% negative mu part
Lnegmumat = 1/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
Lneggmat = diag(ones(nx,1)) + diag(ones(nx-1,1),1);

%% H1
G = cell(1,2);
G{1} = reshape(Lposmumat,[1,size(Lposmumat),1]);
mupos = diag([mu(1:nmu/2), zeros(1,nmu/2)]);
G{2} = reshape(mupos,[1,nmu,nmu,1]);
ttLposmu = cell2core(tt_matrix, G);
%convert to qtt
qttLposmu = convert_H_to_qtt_1_group_fn(ttLposmu,tol);

% Lposg
G{1} = reshape(Lposgmat,[1,size(Lposgmat),1]);
G{2} = G{2}~=0;
ttLposg = sigma/2*cell2core(tt_matrix, G);
%convert to qtt
qttLposg = convert_H_to_qtt_1_group_fn(ttLposg,tol);

qttH1 = qttLposmu + qttLposg;


%% H2

G = cell(1,2);
G{1} = reshape(Lnegmumat,[1,size(Lnegmumat),1]);
muneg = diag([zeros(1,nmu/2), mu(nmu/2+1:nmu)]);
G{2} = reshape(muneg,[1,nmu,nmu,1]);
ttLnegmu = cell2core(tt_matrix, G);
qttLnegmu = convert_H_to_qtt_1_group_fn(ttLnegmu,tol);

% only update G{1} and G{3}
G{1} = reshape(Lneggmat,[1,size(Lneggmat),1]);
G{2} = G{2}~=0;
ttLnegg = sigma/2*cell2core(tt_matrix, G);
qttLnegg = convert_H_to_qtt_1_group_fn(ttLnegg,tol);

qttH2 =qttLnegmu + qttLnegg;

%forming H in qtt
qttH = qttH1 + qttH2;

%% Construct the scattering operator Hs

% Hs positive
tempXpos = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1));
tempXpos_noBC = tempXpos;
tempXpos_noBC(1,:) = 0;

tempmupos = [ones(nmu/2,1);zeros(nmu/2,1)]*wgt';

G = cell(1,2);
G{1} = reshape(tempXpos_noBC,[1,nx,nx,1]);
G{2} = reshape(tempmupos,[1,nmu,nmu,1]);
ttHspos = sigma_s*cell2core(tt_matrix,G);
qttHspos = convert_H_to_qtt_1_group_fn(ttHspos,tol);

% Hs negative
tempXneg = permute(tempXpos,[2,1]);
tempXneg_noBC = tempXneg;
tempXneg_noBC(end,:) = 0;

G = cell(1,2);
G{1} = reshape(tempXneg_noBC,[1,nx,nx,1]);
G{2} = reshape(flipud(tempmupos),[1,nmu,nmu,1]);
ttHsneg = sigma_s*cell2core(tt_matrix,G);
qttHsneg = convert_H_to_qtt_1_group_fn(ttHsneg,tol);

% ttHs = ttHspos + ttHsneg;
qttHs = qttHspos + qttHsneg;

%% Construct the fission Hf
qttHf = qttHs*(chi*nusigma_f)/sigma_s;

%% Fixed-point scheme in tensor train
% fixed point schemes
nxqtt = 2*ones(1,log2(nx));
nmuqtt = 2*ones(1,log2(nmu));

qttvecsize = [nxqtt,nmuqtt]';

qttPsi0 = tt_ones(qttvecsize,numel(qttvecsize),1);
k = 1;
niter = param.niter;
convrate = nan(1,niter);
Itertime = zeros(1,niter);

for iter = 1:niter
  tic
  % fprintf('Iteration = %d/%d \n convrate = %.5e \n elapsed time = %.2f s \n', ...
    % iter,niter, convrate(max(iter-1,1)), sum(Itertime) );
  RHS = amen_mv(qttHs,qttPsi0,tol,'verb',0) + 1/k*amen_mv(qttHf,qttPsi0,tol,'verb',0);

  %  Solve a linear system for Psi1 = Htt\RHS;
  amen_solve_opts = {'verb',0,'ismex',0,'resid_damp',0.1,'x0',0};
  amen_solve_opts{end} = round(qttPsi0, tol);
  qttPsi1 = amen_solve2(qttH,RHS,tol, amen_solve_opts);

  %update eigenvalue
  k = k*sum(amen_mv(qttHf,qttPsi1,tol,'verb',0))./...
    sum(amen_mv(qttHf,qttPsi0,tol,'verb',0));

  %relative error
  convrate(iter) = norm(qttPsi1-qttPsi0)/norm(qttPsi0);
  conv_test2 = 0;
  if iter>5
    conv_test2 = convrate(iter)>convrate(iter-5);
  end
  if (convrate(iter) < fixed_point_tol) || conv_test2
    fprintf('The fixed point scheme converged at iter = %d \n', iter)
    break
  end
  qttPsi0 = qttPsi1;
  Itertime(iter) = toc;
end
%normalize ttPsi1
qttPsi1 = qttPsi1/norm(qttPsi1);
end

