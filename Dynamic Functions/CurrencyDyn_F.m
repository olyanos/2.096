function [dxdt] = CurrencyDyn_F(x,Params,dWdt)
% CurrencyDyn_F   Compute derivatives

% Inputs:         1) x:      state inputs
%                 2) Params: system parameters
%                 3) dWdt:   stochastic input
%
% Outputs:        1) dxdt:   state derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse state vector
n  = Params.n;
Y  = x(1:n);
Yt = x((n+1):2*n);
mu = x((2*n+1):3*n);

% Parse "Params" vector
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

% Dynamic parameters
gam_damp = Params.gam_damp;
alpha    = Params.alpha;
sigma    = Params.sigma;
tau_1    = Params.tau_1;
tau_2    = Params.tau_2;
Yref     = Params.Yref;

% Compute BTRs
rho  = 1;
[Rx] = Build_BLT_Rx(gam1,gam2,d,nu,Yt,rho,E);

% Compute Flows
[x] = Build_x_EXP(y,yw,P,Rx,Gamma,Lambda);
[q] = A*x;
[N] = Lambda*x - (Lambda*((Gamma*Yt).*q))./Yt;

% Compute Derivatives
dY_dt  = (mu + sigma.*dWdt).*Y + gam_damp.*(Yref - Y);
dYt_dt = (Y-Yt)./tau_1;
dmu_dt = (alpha.*N - mu)./tau_2;

% Full vector
dxdt = [dY_dt; dYt_dt; dmu_dt];

end

