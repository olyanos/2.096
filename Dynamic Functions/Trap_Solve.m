function [xm] = Trap_Solve(Fp,Jp,t_vec,dt,x0,u)

% Initialize matrix for results
xm      = zeros(length(x0), length(t_vec));
xm(:,1) = x0;

% Loop over time vector
for jj = 2:length(t_vec)
    
    % Solve Newton
    x0   = xm(:,jj-1);
    u0   = u(:,jj-1);
    u1   = u(:,jj);
    [x1] = Newton_Trap(x0,dt,Fp,Jp,u0,u1);
    
    % Increment
    xm(:,jj) = x1;
end

end

% Newton
function [x1] = Newton_Trap(x0,dt,Fp,Jp,u0,u1)
tol      = 1e-8;
maxIters = 500;
x1       = x0;

% Newton loop
for iter = 1:maxIters
    
    % Define function FT
    F0 = Fp(x0,u0);
    F1 = Fp(x1,u1);
    FT = x1 - x0 - (dt/2)*(F0 + F1);
    
    % Define Jacobian JT
    J  = Jp(x1,u1);
    JT = eye(length(x1)) - (dt/2)*J;
    
    % Solve system
    dx    = -JT\FT;
    nf    = norm(FT);            % norm of f at step k+1
    x1    = x1+dx;               % solution x at step k+1

    % Test
    if nf  < tol           % check for convergence
        % check for convergence
        % fprintf('Converged in %d iterations\n',iter);
        break; 
    end
end

if iter==maxIters % check for non-convergence
    % fprintf('Non-Convergence after %d iterations!!!\n',iter); 
end
end