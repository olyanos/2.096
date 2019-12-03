%% Perform TD Simulation
clear; clc;

%% 1. Run steady-state solver
Run_MLSS
nu  = Params.nu;  % Save these
GDP = Params.y;   % Save these

%% 2. Define dynamic system parameters
Params.gam_damp = 0;                 % Interest damping - ignore
Params.alpha    = 0.1*ones(n,1);     % Converts imbalance to drift
Params.sigma    = 1;                 % Volatility parameter
Params.tau_1    = 1;                 % 1 day is the response lag btwn Y and Yt
Params.tau_2    = 2;                 % 2 day response lag btwn N and mu
Params.Yref     = 0;                 % Ignore for now
dt              = 1/24;              % Time step
tf              = 24*5;              % Run time
t_vec           = 0:dt:tf;           % Time vector

%% 3. Initialize steady state (mu = 0, Y = Yt = Y_steady_state)
mu = zeros(n,1);
x0 = [Yt; Yt; mu];
h  = 1e-5;

%% 4. Anonymize Function and Jacobian
[F_CD] = @(x,dWdt)CurrencyDyn_F(x,Params,dWdt);  % Currency Dynamics Function
[J_CD] = @(x,dWdt)CurrencyDyn_J(F_CD,x,dWdt,h);  % Currency Dynamics Jacobian

%% 5. Perform Time Domain Integration #1: Steady State
%     This is a steady state simulation: there are no inputs,
%     so variables shouldn't move at all.
dWdt_vec = zeros(n,length(0:dt:tf));   % Define a stochastic input
[xm1]    = Trap_Solve(F_CD,J_CD,t_vec,dt,x0,dWdt_vec);

%% 6. Perform Time Domain Integration #2: "Venezuela"
%     In this simulation, a tenth of the way through, import resistance
%     increases but export resistance drops
tf                 = 365;
t_vec1             = 0:dt:round(tf/10);
t_vec2             = (round(tf/10)+dt):dt:tf;

% No stochastic inputs
dWdt_vec1          = zeros(n,length(t_vec1));
dWdt_vec2          = zeros(n,length(t_vec2));
Params.nu          = nu;
Params.y           = GDP;
[F_CD]             = @(x,dWdt)CurrencyDyn_F(x,Params,dWdt);
[J_CD]             = @(x,dWdt)CurrencyDyn_J(F_CD,x,dWdt,h);
[xm2a]             = Trap_Solve(F_CD,J_CD,t_vec1,dt,x0,dWdt_vec1);

% Alter "nu" paramerers
v1                 = 10:9:82;
Params.nu(1:9)     = 1.5*nu(1:9);  % Export resistance
Params.nu(v1)      = 0.5*nu(v1);   % Import resistance
[F_CD]             = @(x,dWdt)CurrencyDyn_F(x,Params,dWdt);  % Redefine
[J_CD]             = @(x,dWdt)CurrencyDyn_J(F_CD,x,dWdt,h);  % Redefine
[xm2b]             = Trap_Solve(F_CD,J_CD,t_vec2,dt,x0,dWdt_vec2);

%% Parse
r.Y  = [xm2a(1 :10,:) xm2b(1 :10,:)];
r.Yt = [xm2a(11:20,:) xm2b(11:20,:)];
r.mu = [xm2a(21:30,:) xm2b(21:30,:)];
r.tv = [t_vec1 t_vec2];

% Initialize
r.N_mat = zeros(10,length(r.tv));
r.q_mat = zeros(90,length(r.tv));
r.x_mat = zeros(90,length(r.tv));

% Use Results to compute tradeflows and N (trade imbalance)
for ii = 1:length(r.tv)
    time_step = r.tv(ii);
    
    % What are the value differential parameters
    if time_step < t_vec2(1)
        Params.nu = nu;
    else
        Params.nu(1:9) = 1.5*nu(1:9);  % Export resistance
        Params.nu(v1)  = 0.5*nu(v1);   % Import resistance
    end
    
    % Compute flows
    Ytv    = r.Yt(:,ii);
    E      = Params.E;
    A      = Params.A;
    Lambda = Params.Lambda;
    Gamma  = Params.Gamma;
    y      = Params.y;
    P      = Params.P;
    d      = Params.d;
    yw     = Params.yw;
    gam1   = Params.gam1;
    gam2   = Params.gam2;
    
    % Compute Flows
    [Rx] = Build_BLT_Rx(gam1,gam2,d,Params.nu,Ytv,rho,E);
    [x]  = Build_x_EXP(y,yw,P,Rx,Gamma,Lambda);
    [q]  = A*x;
    [N]  = Lambda*x - (Lambda*((Gamma*Ytv).*q))./Ytv;
    
    % Save data
    r.x_mat(:,ii) = x;
    r.q_mat(:,ii) = q;
    r.N_mat(:,ii) = N;
end

% Save data in the structure "r"
save('Venezuela_sim_data','r')

%% Plot Currencies
plot(0:dt:tf,[xm2a(1:10,:) xm2b(1:10,:)])
set(gca,'FontName','Times','FontSize',13)
xlabel('Time (hours)')
ylabel('Currency ($)')
legend('Australia', 'Brazil', 'Canada', 'EU', 'Ghana', 'Israel', 'Malysia', 'Singapore', 'Switzerland', 'Turkey')
set(gcf,'Units','inches','Position',[0 0 9 3.5])
tightfig(gcf) % Now "File" => "Export Setup" => "Expand axes to fill figure"

%% 7. Perform Time Domain Integration #3: Stochastic Input
tf                 = 30;                % Run time
t_vec              = 0:dt:tf;           % Time vector
Params.nu          = nu;
Params.y           = GDP;
[F_CD]             = @(x,dWdt)CurrencyDyn_F(x,Params,dWdt);
[J_CD]             = @(x,dWdt)CurrencyDyn_J(F_CD,x,dWdt,h);
dWdt_vec           = [zeros(10,1) 0.01*randn(n,length(t_vec)-1)];
[xm3]              = Trap_Solve(F_CD,J_CD,t_vec,dt,x0,dWdt_vec);
