clear all; clc;
tic
% Nodes
n = 10;
m = n-1;
% Australia, Canada, China, EU, India, Chile, Japan, Mexico, Morocco, Poland
[Y0_real] = [0.68;0.75;0.14;1.11;0.14;0.13;0.92;0.52;0.10;0.26];
% Define system parameters and add them to "Params" structure
Params.n    = n;             % number of nodes
Params.m    = n-1;           % number of export flows on each node
Params.y    = [1.43; 1.71; 13.6; 18.8; 2.73;0.298;4.97;1.22;0.118;0.586];     % nodal GDPs (billions)
Params.P    = [117.898;114.525;121.559;111.247;167.598;128.6239518;...
               104.9814383;136.5766487;110.8502192;111.6253484];     % nodal CPIs

fileID = fopen('distvec.txt', 'r');
formatSpec = '%f';
dist = fscanf(fileID, formatSpec);
Params.d    = reshape(dist,10,9);   % distance vector (km)          

Params.yw   = sum(Params.y);          % global GDP

Params.gam1 = ones(n*m,1);   % tuning vector 1
Params.gam2 = ones(n*m,1).*0.25;   % tuning vector 2
Array_errors = ones(1,1)*100;
Array_errors2 = ones(1,1)*100;
error2= 1e6;
error = 1e6;
while(error2> 0.5 | isnan(error2) )
    while (error > 0.6)
        Params.nu = rand(10,9);
        % Build special incidence matrix E, Lambda, Gamma and A
        [E,Lambda,Gamma] = Build_ELG(Params.n);
        [A]              = Build_A(Params.n);
        % Solve linear system
        [Y0] = MLSS_Lin_Solve(Gamma,Lambda,A,Params);
        lin_solv_diff_from_real = abs(Y0 - Y0_real)./Y0_real;
        error = norm(lin_solv_diff_from_real,inf);
        Array_errors = [Array_errors;error];
    end
    %% Use loop to drive N to 0 for a given rho:
    Yt  = Y0;
    h = 1e-5; % finite diff paramter for numerical Jacobian.
    for rho = 0.1 : 0.1 : 1     
        [N] = MLSS_NonLin(Yt,Gamma,Lambda,A,Params,rho);
        [newton_sol,k] = Newton(Yt,Gamma,Lambda,A,Params,rho,h);        
        Yt = newton_sol;
    end
    diff_from_origianl = abs(Y0_real - Yt)./Y0_real;
    error;
    error2= norm(diff_from_origianl,inf);
    Array_errors2 = [Array_errors2;error2];
    error = 1e6;
end
toc