%% Solve for Steady State Currency Values
clear variables; clc

% Nodes
n = 5;
m = n-1;
% Australia, Canada, China, EU, India
% Define system parameters and add them to "Params" structure
Params.n    = n;             % number of nodes
Params.m    = n-1;           % number of export flows on each node
Params.y    = [1.43; 1.71; 13.6; 18.8; 2.73];     % nodal GDPs (billions)
Params.P    = [118; 115; 122; 111; 168];     % nodal CPIs
Params.d    = [16000; 9000; 16700; 10300;
                16000; 10000; 6000; 11000;
                9000; 10000; 8000; 4000;
                17000; 6000; 8000; 6000;
                10000; 11000; 4000; 6000];   % distance vector
Params.nu   = [0.0482; 0.2327; 0.3051; 0.4518; 
                1.0006; 0.8715; 0.1926; 1.9027;
                1.2965; 1.8060; 0.8578; 0.3376;
                1.5178; 0.6343; 0.4766; 0.7606;
                1.1529; 0.1558; 1.5679; 0.3409];   % value differential vector
Params.yw   = 84.4;          % global GDP
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
[newton_sol,k] = Newton(Yt,Gamma,Lambda,A,Params,rho,h);
fprintf('number of Newton iterations=%d\n',k);
