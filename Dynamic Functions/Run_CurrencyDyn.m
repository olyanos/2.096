%% Time Domain Integration
clear; clc

% System Definition
n = 5;
m = n-1;

% Define system parameters and add them to "Params" structure
Params.n    = n;               % number of nodes
Params.m    = n-1;             % number of export flows on each node
Params.y    = rand(n,1);       % nodal GDPs
Params.P    = rand(n,1);       % nodal CPIs
Params.d    = rand(n*m,1);     % distance vector
Params.nu   = rand(n*m,1);     % value differential vector
Params.yw   = rand;            % global GDP
Params.gam1 = 0.1*ones(n*m,1); % tuning vector 1
Params.gam2 = ones(n*m,1);     % tuning vector 2

% Build special incidence matrix E, Lambda, Gamma and A
[E,Lambda,Gamma] = Build_ELG(Params.n);
[A]              = Build_A(Params.n);

% Load into param structure
Params.E      = E;
Params.A      = A;
Params.Lambda = Lambda;
Params.Gamma  = Gamma;

% Dynamic parameters
Params.gam_damp = 0*ones(n,1);       % Interest damping - ignore
Params.alpha    = 0.1*ones(n,1);     % Converts imbalance to drift
Params.sigma    = 0.1;               % Volatility parameter
Params.tau_1    = 0.2;                % 12 hours is the response lag btwn Y and Yt
Params.tau_2    = 0.2;              % 1 week response lag btwn N and mu
Params.Yref     = 0;                 % Ignore for now

%%%% --- need Omar's nonlinear solution here --- %%%%
Y  = rand(n,1); % This should be the stead state solution
Yt = Y;
mu = zeros(n,1);
x0 = [Y; Yt; mu];
h  = 1e-5;

% Anonymize Function and Jacobian
[F_CD] = @(x,dWdt)CurrencyDyn_F(x,Params,dWdt);  % Currency Dynamics Function
[J_CD] = @(x,dWdt)CurrencyDyn_J(F_CD,x,dWdt,h);  % Currency Dynamics Jacobian

% Define a stochastic input
dt       = 0.1;     % time step (1 hour)
tf       = 1;      % final time (5 days)

% Define a stochastic input - Ignore all of this
OU_tau   = 1e7;   % stochastic param
OU_alpha = 1e3;   % stochastic param
nLu      = 1e5;   % stochastic param
dWdt     = zeros(n,length(0:dt:tf));
for ii = 1:n
    [OU_ut,OU_Time] = OU_Sim(dt,tf,OU_tau,OU_alpha);
    dWdt(ii,:)      = 1e5* OU_ut;
end

% Time parameterization
% puT   = spline(OU_Time,1e5*nV);
% dW_dt = @(t)ppval(puT,t);

% Perform Time Domain Integration
t_vec = 0:dt:tf;
[xm]  = Trap_Solve(F_CD,J_CD,t_vec,dt,x0,dWdt);



