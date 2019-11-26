function [Y0] = MLSS_Lin_Solve(Params)
%BUILD_MLSS       Build the steady state Marshall-Lerner function vector

% Inputs:         1) Params: system structure and parameters
%
% Outputs:        1) Y0:     linear system solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0) Parse "Params" vector
E      = Params.E;
A      = Params.A;
Lambda = Params.Lambda;
Gamma  = Params.Gamma;
n      = Params.n;
m      = Params.m;
y      = Params.y;
P      = Params.P;
d      = Params.d;
nu     = Params.nu;
yw     = Params.yw;
gam1   = Params.gam1;
gam2   = Params.gam2;

% 1) Set rho = 0;
rho = 0;

% 2) Get BTRs for rho = 0
Yt   = zeros(n,1); % This is arbitrary
[Rx] = Build_BLT_Rx(gam1,gam2,d,nu,Yt,rho,E);

% 3) Compute static flows
[x]  = Build_x_EXP(y,yw,P,Rx,Gamma,Lambda);

% 4) Build and partition matrix M
M  = diag(Lambda*x) - Lambda*diag(A*x)*Gamma;
M1 = M(1,1);
M2 = M(1,2:end);
M3 = M(2:end,1);
M4 = M(2:end,2:end);

% 5) Assign and solve
x1 = 1.30;
xn = -M4\M3*x1;
Y0 = [x1; xn];

end

