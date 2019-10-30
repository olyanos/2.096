%% Solve for Steady State Currency Values
clear variables; clc

% Nodes
n = 4;
m = n-1;

% Define system parameters and add them to "Params" structure
Params.n    = n;             % number of nodes
Params.m    = n-1;           % number of export flows on each node
Params.y    = rand(n,1);     % nodal GDPs
Params.P    = rand(n,1);     % nodal CPIs
Params.d    = rand(n*m,1);   % distance vector
Params.nu   = rand(n*m,1);   % value differential vector
Params.yw   = rand;          % global GDP
Params.gam1 = ones(n*m,1);   % tuning vector 1
Params.gam2 = ones(n*m,1);   % tuning vector 2

% Build special incidence matrix E, Lambda, Gamma and A
[E,Lambda,Gamma] = Build_ELG(Params.n);
[A]              = Build_A(Params.n);

%% Solve linear system
[Y0] = MLSS_Lin_Solve(Gamma,Lambda,A,Params);

%% Use loop to drive N to 0 for a given rho:
rho = 0.1; % scaling newton
Yt  = Y0;
[N] = MLSS_NonLin(Yt,Gamma,Lambda,A,Params,rho);

h = 1e-2;
[newton_sol,k] = Newton(Yt,Gamma,Lambda,A,Params,rho,h)
fprintf('number of Newton iterations=%d\n',k);
