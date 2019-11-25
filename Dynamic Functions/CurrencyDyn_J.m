function [J] = CurrencyDyn_J(F,x,dWdt,h)
% CurrencyDyn_J   Compute derivatives

% Inputs:         1) F:      NL function
%                 2) x:      state inputs
%                 3) dWdt:   stochastic input
%                 4) h:      perturbation size
%
% Outputs:        1) J:      numerical Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n  = length(x);
J  = zeros(n,n);
Y0 = F(x,dWdt);
for k = 1:n
    x_new   = x + h;
    Y0_new  = F(x_new,dWdt);
    J(:,k)  = (Y0_new-Y0)./h;
end
end

