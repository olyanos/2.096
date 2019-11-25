function [OU_Noise,OU_Time] = OU_Sim(dt,tf,OU_tau,OU_alpha)
% OU_SIM          Generate an OU vector of noise
%
% Inputs:         
% 1) dt           Time step of simulation
% 2) tf           Stopping time
% 3) OU_tau       OU time constant
% 4) OU_alpha     OU scalar
%
% u_dot = (-OU_alpha*u + randn)/OU_tau
%
% Outputs:
% 1) OU_Noise     OU noise vector of length 0:dt:tf
% 2) OU_Time      OU noise time vector
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OU Noise Variable
syms u(t) eta(t)

% DAE Variables
ODEvars = u(t);

% DAEs
ODEs    = diff(u(t))==(-OU_alpha*u(t) + eta(t))/OU_tau;

% Set up system  
[Mm,Ff] = massMatrixForm(ODEs,ODEvars);
Mm      = odeFunction(Mm,ODEvars);
Ff      = odeFunction(Ff,ODEvars, eta(t));

% Build a vector of Random Noise Values
N_vec = [0; randn(length(0:dt:tf)-1,1)];

% Use a spline add noise
t    = linspace(0,tf,length(N_vec));
pEta = spline(t,N_vec);

% Voltage Magnitude and Phase Ocillation
eta = @(t) ppval(pEta,t);
Ff  = @(t,Y) Ff(t,Y,eta(t));

% Set options
opt = odeset('mass',Mm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan         = 0:dt:tf;
[t_out,y_out] = ode45(Ff,tspan,0,opt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! Parse Output !!! %
OU_Time  = t_out;
OU_Noise = y_out;

end
