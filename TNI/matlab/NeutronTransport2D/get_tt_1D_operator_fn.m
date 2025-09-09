function [ ttH] = get_tt_1D_operator_fn (sigma, a, b, nmu, nx)

nE = numel(sigma); % index g
% nmu = 2; % index l
% nx = 3; % index i - number of EDGES
dx = (b-a)/(nx-1); %spacial step

% regenerate mu and weight
[mu,wgt] = legpts(nmu,[-1,1]);
mu = flip(mu);
wgt = flip(wgt'/2);

%% forming the system in the tensor format
% positive mu part
Lposmumat = mu(1)/dx*(diag(ones(nx,1)) - diag(ones(nx-1,1),-1));
Lposgmat = diag(ones(nx,1)) + diag(ones(nx-1,1),-1);

% negative mu part
Lnegmumat = mu(2)/dx*(-diag(ones(nx,1)) + diag(ones(nx-1,1),1));
Lneggmat = diag(ones(nx,1)) + diag(ones(nx-1,1),1);

%% H1

G = cell(1,3);
G{1} = reshape(Lposmumat,[1,size(Lposmumat),1]);
G{2} = reshape([1,0;0,0],[1,2,2,1]);
G{3} = reshape(eye(nE),[1,size(eye(nE)),1]);
ttLposmu = cell2core(tt_matrix, G);

% only update G{1} and G{3}
G{1} = reshape(Lposgmat,[1,size(Lposgmat),1]);
G{3} = reshape(diag(sigma/2),[1,nE,nE,1]);
ttLposg = cell2core(tt_matrix, G);
ttH1 = ttLposmu + ttLposg;

%% H2

G = cell(1,3);
G{1} = reshape(Lnegmumat,[1,size(Lnegmumat),1]);
G{2} = reshape([0,0;0,1],[1,2,2,1]);
G{3} = reshape(eye(nE),[1,size(eye(nE)),1]);
ttLnegmu = cell2core(tt_matrix, G);

% only update G{1} and G{3}
G{1} = reshape(Lneggmat,[1,size(Lneggmat),1]);
G{3} = reshape(diag(sigma/2),[1,nE,nE,1]);
ttLnegg = cell2core(tt_matrix, G);

ttH2 = ttLnegmu + ttLnegg;

%forming H in tt
ttH = ttH1 + ttH2;
