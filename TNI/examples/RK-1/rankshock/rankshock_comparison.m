function output = rankshock_comparison(mname, dtvals, N, reltol, abstol)

addpath(genpath('../../../matlab/RK-1/'))
addpath(genpath('../../../matlab/utils/tt-toolbox/'))
addpath(genpath('../../../matlab/utils/ttfunc/'))

% RK methods (order/embedded order - extrapolation used)
% mname = 'Bogacki-Shampine-ERK'; % 3/2 method
% mname = 'ERK-3-3'; % 3/2 method
% mname = 'Fehlberg-ERK'; % 5/4
% mname = 'Verner-6-5-ERK'; % 6/5

a = 1; b = -3; % A matrix entries
A = diag(b*ones(1,N)) + diag(a*ones(1,N-1),1) + diag(a*ones(1,N-1),-1);

phi = @(j) cos(2*pi*(1:N)'*j/N);
psi = @(j) sin(2*pi*(1:N)'*j/N);
sig = @(j) (3/4).^j;
rlow = 6;        % low rank
rhigh = 25;      % high rank
vlow = zeros(N); vhigh=vlow;  % initialize forcing matrices

%% Forcing matrices for full grid
for j = 1:rlow
    vlow = vlow + phi(j)*psi(j)';
end

for j = 1:rhigh
    vhigh = vhigh + sig(j)*psi(j)*phi(j)';
end
vt = @(t) vlow*(t<5 || t>15) + vhigh*(t>=5 && t<=15); % forcing function
V = {vlow,vhigh};

%% Forcing matrices for TT
G{1} = reshape(phi(1:rlow),1,N,rlow);
G{2} = reshape(psi(1:rlow)',rlow,N,1);
Vlowtt = cell2core(tt_tensor,G);           % lower rank

G{1} = reshape(sig(1:rhigh).*psi(1:rhigh),1,N,rhigh);
G{2} = reshape(phi(1:rhigh)',rhigh,N,1);
Vhightt = cell2core(tt_tensor,G);           % higher rank
Vtt = {Vlowtt,Vhightt};

%% Right-hand side functions and initial conditions
frhs = @(t,Ytt,et) tt_oderhs(t,Ytt,a,b,Vtt,et);  % tt right hand side
oderhs = @(t,F) fg_oderhs(t,F,a,b,V);

F0 = zeros(size(A));    % initial condition
F0tt  = tt_zeros(N,2);

%% Time, Tolerances, RK method
tvals  = 0:0.5:20;      % time samples for solution comparison

dtinit = dtvals(1); dtmin = dtvals(2); dtmax = dtvals(3);

rtol = reltol;               % relative tolerance for tt RK solver
atol = abstol*tt_ones(N,2); % absolute tolerance for tt/fg RK solver
atolfg = abstol*ones(N);    % absolute tolerance for fg RK solver
ett = rtol;                % tt-rank truncation tolerance

B = butcher(mname);
s = numel(B(1,:))-1;

%% Compute reference solution
fprintf('\nCalculating a reference solution ...\n')
tstart = tic;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[tout,Fapprox] = ode45(@(t,F) full_oderhs(t,F,A,V),tvals,F0(:),opts);
trefsol = toc(tstart);
fprintf('Thanks for your patience. That took %g seconds. \n',trefsol)

% Storage initialization
refrank = zeros(size(tout));
Fttref = cell(size(tout));
FE = cell(size(tout));

% Get reference ranks
for i = 1:length(tout)
    Fttref{i} = tt_tensor(Fapprox(i,:),ett,[N,N]);
    refrank(i) = max(Fttref{i}.r);
    FE{i} = reshape(Fapprox(i,:),N,N);
end

%% Full grid ERK solve
fprintf('\nRunning with FG-ERK integrator: %s (order = %i)\n',mname,B(s+1,1))

par.B = B;
par.rtol = rtol;
par.dtmin = dtmin;
par.dtmax = dtmax;
par.dtinit = dtinit;
par.atol = atolfg;

tstart = tic;
[~,FA,nstepsRK,fgoutput] = solve_ERK_struct(oderhs,tvals,F0,par);
fg_time = toc(tstart);
fprintf('FG time elapsed: %g s\n',fg_time)
fprintf('FG nsteps:       %g  \n',nstepsRK)
fprintf('FG failed steps: %g  \n',fgoutput.nfail)
fprintf('FG Failure rate: %.2f%% \n', fgoutput.nfail/nstepsRK*100)


% Storage initialization
fgerror = zeros(size(tout));
fgrmserror = zeros(size(tout));

for i = 1:length(tout)
    error  = FE{i}-FA{i};
    fgerror(i) = norm(error,'fro')/norm(FE{i},'fro');
    fgrmserror(i) = norm(error.*error/N^2,'fro');
end

%% TT ERK solve
par.B = B;
par.rtol = rtol;
par.dtmin = dtmin;
par.dtmax = dtmax;
par.dtinit = dtinit;
par.atol = atol;
par.ettinit = ett;

fprintf('\nRunning with TT-ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
tstart = tic;
[~,FttRK23,nstepstt,ttoutput] = solve_ERKtt_struct(frhs,tvals,F0tt,par);
tt_time = toc(tstart);
fprintf('TT time elapsed: %g s \n',tt_time)
fprintf('TT nsteps:       %g  \n',nstepstt)
fprintf('TT failed steps: %g  \n',ttoutput.nfail)
fprintf('TT Failure rate: %.2f%% \n', ttoutput.nfail/nstepstt*100)
fprintf('-------------------------------------------\n')

% Storage initialization
tterror = zeros(size(tout));
maxranks = tterror;

% Get error and ranks
for i = 1:length(tout)
    error  = Fttref{i} - FttRK23{i};
    tterror(i) = norm(error)/norm(Fttref{i});         % relative errror
    maxranks(i) = max(FttRK23{i}.r);
end

output.mname = mname;
output.N = N;
output.rtol = rtol; 
output.abstol = abstol;
output.dtvals = dtvals;
output.ett = ett;
output.tt_time = tt_time;
output.fg_time = fg_time;
output.nstepstt = nstepstt;
output.nstepsfg = nstepsRK;
output.nfailtt = ttoutput.nfail;
output.nfailfg = fgoutput.nfail;
output.ttoutput = ttoutput;
output.ttranks = maxranks;
output.refranks = refrank;
output.tterror = tterror;
output.fgerror = fgerror;

end