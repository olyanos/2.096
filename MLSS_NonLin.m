function [N] = MLSS_NonLin(Yt,Gamma,Lambda,A,Params,rho)
%BUILD_MLSS       Build the steady state Marshall-Lerner function vector

% Inputs:         1) Yt:     effective nodal currency vector
%                 2) Gamma:  selection matrix
%                 3) Lambda: summation matrix
%                 4) A:      matrix mapping q=A*x
%                 5) Params: system parameters
%                 6) rho:    parametrives the nonlinearity
%
% Outputs:        1) N:  trade imbalance (scaled by currency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % --- For a given currency vector Yt --- % % %

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

% 1) Compute BTR
[Rx] = Build_BLT_Rx(gam1,gam2,d,nu,Yt,rho,n,m);

% 2) Compute flows;
[x] = Build_x_EXP(y,yw,P,Rx,n,m);

% 3) Compute N
N = (diag(Lambda*x) - Lambda*diag(A*x)*Gamma)*Yt;

end

