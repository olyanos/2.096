function [N] = MLSS_NonLin(Yt,Params,rho)
%BUILD_MLSS       Build the steady state Marshall-Lerner function vector

% Inputs:         1) Yt:     effective nodal currency vector
%                 2) Params: system parameters
%                 3) rho:    parametrives the nonlinearity
%
% Outputs:        1) N:  trade imbalance (scaled by currency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % --- For a given currency vector Yt --- % % %

% 0) Parse "Params" vector
E      = Params.E;
A      = Params.A;
Lambda = Params.Lambda;
Gamma  = Params.Gamma;
y      = Params.y;
P      = Params.P;
d      = Params.d;
nu     = Params.nu;
yw     = Params.yw;
gam1   = Params.gam1;
gam2   = Params.gam2;

% 1) Compute BTR
[Rx] = Build_BLT_Rx(gam1,gam2,d,nu,Yt,rho,E);

% 2) Compute flows;
[x] = Build_x_EXP(y,yw,P,Rx,Gamma,Lambda);

% 3) Compute N
N = (diag(Lambda*x) - Lambda*diag(A*x)*Gamma)*Yt;

end

