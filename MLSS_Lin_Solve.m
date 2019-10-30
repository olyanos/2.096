function [Y0] = MLSS_Lin_Solve(Gamma,Lambda,A,Params)
%BUILD_MLSS       Build the steady state Marshall-Lerner function vector

% Inputs:         1) Yt:     effective nodal currency vector
%                 2) Gamma:  selection matrix
%                 3) Lambda: summation matrix
%                 4) A:      matrix mapping q=A*x
%                 5) Params: system parameters
%                 6) rho:    parametrives the nonlinearity
%
% Outputs:        1) Y0:     linear system solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0) Parse "Params" vector
n    = Params.n;
m    = Params.m;
y    = Params.y;
P    = Params.P;
d    = Params.d;
nu   = Params.nu;
yw   = Params.yw;
gam1 = Params.gam1;
gam2 = Params.gam2;

% 1) Set rho = 0;
rho = 0;

% 2) Get BTRs for rho = 0
Yt   = zeros(n,1); % This is arbitrary
[Rx] = Build_BLT_Rx(gam1,gam2,d,nu,Yt,rho,n,m);

% 3) Compute static flows
[x]  = Build_x_EXP(y,yw,P,Rx,n,m);

% 4) Build and partition matrix M
M  = diag(Lambda*x) - Lambda*diag(A*x)*Gamma;
M1 = M(1,1);
M2 = M(1,2:end);
M3 = M(2:end,1);
M4 = M(2:end,2:end);

% 5) Assign and solve
x1 = 1;
xn = -M4\M3*x1;
Y0 = [x1; xn];

end

