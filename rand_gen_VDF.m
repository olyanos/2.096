clear all; clc;
tic
% Nodes
n = 10;
m = n-1;
% Australia, Brazil, Canada, EU, Ghana, Israel,Malysia,Singapore,Switzerland,
% Turkey
[Y0_real] = [1.30;3.19;1.30;0.89;4.35;3.60;4.30;1.38;0.98;3.64];
% Define system parameters and add them to "Params" structure
Params.n    = n;             % number of nodes
Params.m    = n-1;           % number of export flows on each node
Params.y    = [1.33;2.05;1.65;17.3;0.06;0.35;0.315;0.338;0.68;0.85];     % nodal GDPs (billions)
Params.P    = [115.6867846;155.6687862;111.9848311;109.3524635;232.2564866;106.3806971;119.6050658;113.2661463;98.26686201;174.9687033];     % nodal CPIs

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
    error2= norm(diff_from_origianl,inf)
    Array_errors2 = [Array_errors2;error2];
    error = 1e6;
end
toc