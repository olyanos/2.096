%% Solve for Steady State Currency Values
clear variables; clc

time = [];
time_jfree = [];
mem = [];
mem_jfree = [];
for i = 1:100
    % Nodes
n = 20;
m = n-1;
% Define system parameters and add them to "Params" structure
 Params.n    = n;             % number of nodes
 Params.m    = n-1;           % number of export flows on each node
 Params.y = rand(n,1);
 Params.P = rand(n,1)*100;
 Params.d = rand(n*m,1)*10000;
 Params.nu = rand(n*m,1);
 Params.yw = rand(1,1)*100;
 
% Australia, Canada, China, EU, India, Chile, Japan, Mexico, Morocco, Poland
% Params.y    = [1.43; 1.71; 13.6; 18.8; 2.73; 
%                 0.298; 4.97; 1.22; 0.118; 0.586];     % nodal GDPs (trillions)
% Params.P    = [117.898; 114.525; 121.559; 111.247; 167.598;
%                 128.624; 104.981; 136.577; 110.850; 111.625];    % nodal CPIs (2010 = 100)
%             
% fileID = fopen('distvec.txt', 'r');
% formatSpec = '%f';
% dist = fscanf(fileID, formatSpec);
% Params.d    = dist;   % distance vector (km)
% 
% Params.nu   = [0.0482; 0.2327; 0.3051; 0.4518; 
%                 1.0006; 0.8715; 0.1926; 1.9027;
%                 1.2965; 1.8060; 0.8578; 0.3376;
%                 1.5178; 0.6343; 0.4766; 0.7606;
%                 1.1529; 0.1558; 1.5679; 0.3409];   % value differential vector
% Params.yw   = sum(Params.y);          % global GDP
% Params.gam1 = ones(n*m,1);   % tuning vector 1
% Params.gam2 = ones(n*m,1);   % tuning vector 2

% Build special incidence matrix E, Lambda, Gamma and A
[E,Lambda,Gamma] = Build_ELG(Params.n);
[A]              = Build_A(Params.n);

%% Solve linear system
[Y0] = MLSS_Lin_Solve(Gamma,Lambda,A,Params);

%% Use loop to drive N to 0 for a given rho:
Yt  = Y0;
h = 1e-5; % finite diff paramter for numerical Jacobian.
rho = 0.05;
    
    [N] = MLSS_NonLin(Yt,Gamma,Lambda,A,Params,rho);
    tic; [newton_sol,k, m] = Newton(Yt,Gamma,Lambda,A,Params,rho,h);
    t = toc;
    fprintf('number of Newton iterations=%d\n',k);

    tic; [newton_sol_jfree,k_jfree, m_jfree] = Newton_JFree(Yt,Gamma,Lambda,A,Params,rho, 1e-02, 1e-5);
    t_jfree = toc;
    fprintf('number of Newton iterations=%d\n',k_jfree);

    time = [time, t];
    time_jfree = [time_jfree, t_jfree];
    mem = [mem, m];
    mem_jfree = [mem_jfree, m_jfree];
    Yt = newton_sol;
end

plot(time)
hold on
plot(time_jfree)
hold off
legend('backslash', 'j free')
ylabel('time')

plot(mem)
hold on
plot(mem_jfree)
hold off
legend('backslash', 'j free')
ylabel('memory')

