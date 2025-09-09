function [Psi1, k, convrate] = fg_1DSlab_EVP_fn(param, tol)

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
%% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]);
mu = flip(mu');
wgt = flip(wgt'/2);

%% forming the system in the tensor format
Lposmu = 1/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));
% Lposmu = diag(ones(nx,1)) - diag(ones(nx-1,1),-1);
mupos = diag([mu(1:nmu/2), zeros(1,nmu/2)]);
Lposmu = ttt(tensor(Lposmu), ttt(tensor(mupos),tensor(eye(nE))));

Lposg = diag(ones(nx,1)) + diag(ones(nx-1,1),-1);

Lposg = ttt(tensor(Lposg), ttt(tensor(mupos~=0), tensor(diag(sigma/2))));
Lpos = Lposmu +  Lposg;

Lnegmu = 1/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
muneg = diag([zeros(1,nmu/2), mu(nmu/2+1:nmu)]);

Lnegmu = ttt(tensor(Lnegmu), ttt(tensor(muneg),tensor(eye(nE))));

Lnegg = diag(ones(nx,1)) + diag(ones(nx-1,1),1);
Lnegg = ttt(tensor(Lnegg), ttt(tensor(muneg~=0), tensor(diag(sigma/2))));
Lneg = Lnegmu + Lnegg;

H = Lpos + Lneg;


%% Construct the scattering operator tensor

tempXpos = tensor(1/2*(diag(ones(nx,1)) + diag(ones(nx-1,1),-1)));
tempXpos_noBC = tempXpos;
tempXpos_noBC(1,:) = 0;
tempmupos = tensor([ones(nmu/2,1);zeros(nmu/2,1)]*wgt');
Hspos = ttt(tempXpos_noBC,tempmupos);

tempXneg = permute(tempXpos,[2,1]);
tempXneg_noBC = tempXneg;
tempXneg_noBC(end,:) = 0;
Hsneg = ttt(tempXneg_noBC, flipud(tempmupos));

Hs = ttt(Hspos + Hsneg, tensor(sigma_s));

%% Construct the fission operator Hf
% resuse Hspos and Hsneg, which is the interpolation and mu_quadrature term
Hf = ttt(Hspos + Hsneg, tensor(chi.*nusigma_f'));

%% Fixed-point scheme in fullgrid
Hsmat = double(reshape(permute(Hs,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]));
Hmat = double(reshape(permute(H,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]));
Hfmat = double(reshape(permute(Hf,[1,3,5,2,4,6]),[nx*nmu*nE,nx*nmu*nE]));
Psi0 = rand(nx*nmu*nE,1);
k = 1;
niter = 100;
convrate = nan(1,niter);
for iter = 1:niter
  RHS = Hsmat*Psi0 + 1/k*Hfmat*Psi0;
  Psi1 = Hmat\RHS;
  k = k*sum(Hfmat*Psi1)/sum(Hfmat*Psi0);
  %   fprintf('||Psi1-Psi0|| = .%5e \n',norm(Psi1-Psi0));
  convrate(iter) = norm(Psi1-Psi0);
  if convrate(iter) < tol
    fprintf('The fixed point scheme converged at iter = %d \n', iter)
    break
  end
  Psi0 = Psi1;
end
%normalize ttPsi1
Psi1 = Psi1/norm(Psi1);



end

