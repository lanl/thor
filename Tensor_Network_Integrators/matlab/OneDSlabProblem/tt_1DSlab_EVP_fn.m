function [ttH, ttHf, ttHs, ttPsi1, k, convrate, Itertime] = tt_1DSlab_EVP_fn(param, tol, savefilename)

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

G = cell(1,3);
G{1} = reshape(Lposmumat,[1,size(Lposmumat),1]);
mupos = diag([mu(1:nmu/2), zeros(1,nmu/2)]);
G{2} = reshape(mupos,[1,nmu,nmu,1]);
G{3} = reshape(eye(nE),[1,size(eye(nE)),1]);
ttLposmu = cell2core(tt_matrix, G);

% Lposg
G{1} = reshape(Lposgmat,[1,size(Lposgmat),1]);
G{2} = G{2}~=0;
G{3} = reshape(diag(sigma/2),[1,nE,nE,1]);
ttLposg = cell2core(tt_matrix, G);
ttH1 = ttLposmu + ttLposg;


%% H2

G = cell(1,3);
G{1} = reshape(Lnegmumat,[1,size(Lnegmumat),1]);
muneg = diag([zeros(1,nmu/2), mu(nmu/2+1:nmu)]);
G{2} = reshape(muneg,[1,nmu,nmu,1]);
G{3} = reshape(eye(nE),[1,size(eye(nE)),1]);
ttLnegmu = cell2core(tt_matrix, G);

% only update G{1} and G{3}
G{1} = reshape(Lneggmat,[1,size(Lneggmat),1]);
G{2} = G{2}~=0;
G{3} = reshape(diag(sigma/2),[1,nE,nE,1]);
ttLnegg = cell2core(tt_matrix, G);

ttH2 = ttLnegmu + ttLnegg;

%forming H in tt
ttH = ttH1 + ttH2;

%% Construct the scattering operator Hs

% Hs positive
tempXpos = 1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1));
tempXpos_noBC = tempXpos;
tempXpos_noBC(1,:) = 0;

tempmupos = [ones(nmu/2,1);zeros(nmu/2,1)]*wgt';

G = cell(1,3);
G{1} = reshape(tempXpos_noBC,[1,nx,nx,1]);
G{2} = reshape(tempmupos,[1,nmu,nmu,1]);
G{3} = reshape(sigma_s,[1,nE,nE,1]);
ttHspos = cell2core(tt_matrix,G);

% Hs negative
tempXneg = permute(tempXpos,[2,1]);
tempXneg_noBC = tempXneg;
tempXneg_noBC(end,:) = 0;

G = cell(1,3);
G{1} = reshape(tempXneg_noBC,[1,nx,nx,1]);
G{2} = reshape(flipud(tempmupos),[1,nmu,nmu,1]);
G{3} = reshape(sigma_s,[1,nE,nE,1]);
ttHsneg = cell2core(tt_matrix,G);
ttHs = ttHspos + ttHsneg;


%% Construct the fission Hf
G = core2cell(ttHspos);
G{3} = reshape(chi.*nusigma_f',[1,nE,nE,1]);
ttHfpos = cell2core(tt_matrix, G);

G = core2cell(ttHsneg);
G{3} = reshape(chi.*nusigma_f',[1,nE,nE,1]);
ttHfneg = cell2core(tt_matrix, G);

ttHf = ttHfpos + ttHfneg;

%% Fixed-point scheme in tensor train
% fixed point schemes
ttPsi0 = tt_ones([nx;nmu;nE],3,1);
k = 1;
niter = param.niter;
convrate = nan(1,niter);
Itertime = zeros(1,niter);

for iter = 1:niter
  tic
  fprintf('Iteration = %d/%d \n convrate = %.5e \n elapsed time = %.2f s \n', ...
    iter,niter, convrate(max(iter-1,1)), sum(Itertime) );
  RHS = amen_mv(ttHs,ttPsi0,tol) + 1/k*amen_mv(ttHf,ttPsi0,tol);

  %  Solve a linear system for Psi1 = Htt\RHS;
  amen_opts = {'x0',round(ttPsi0,tol)};
  ttPsi1 = amen_solve2(ttH,RHS,tol, amen_opts);

  %update eigenvalue
  k = k*sum(amen_mv(ttHf,ttPsi1,tol))./sum(amen_mv(ttHf,ttPsi0,tol));
  %   fprintf('iter = %d - ||Psi1 - Psi0|| = %.5e \n',iter, norm(ttPsi1-ttPsi0));
  %relative error
  convrate(iter) = norm(ttPsi1-ttPsi0)/norm(ttPsi0);
  conv_test2 = 0;
  if iter>5
    conv_test2 = convrate(iter)>convrate(iter-5);
  end
  if (convrate(iter) < fixed_point_tol) || conv_test2
    fprintf('The fixed point scheme converged at iter = %d \n', iter)
    break
  end
  ttPsi0 = ttPsi1;
  Itertime(iter) = toc;

  if nargin>2
    save(savefilename,'ttPsi1','k','convrate',...
  'param','tol','iter');
  end

end
%normalize ttPsi1
ttPsi1 = ttPsi1/norm(ttPsi1);
end

