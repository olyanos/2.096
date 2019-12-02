function [newton_sol,k] = Newton(Yt,Gamma,Lambda,A,Params,rho,h)
% find the 

k = 0;
n = length(Yt);
alpha = 1;
fYt = MLSS_NonLin(Yt,Params,rho);
dYt = ones(n,1);
error1 = 1;
while ( norm (dYt) > 1e-8 & norm (fYt) > 1e-8)
         fYt = MLSS_NonLin(Yt,Params,rho); % function calculated @ (Yt);
         JYt = compute_Jacobian(Yt,Gamma,Lambda,A,Params,rho,h); % Jacobian calculated @ (Yt);
             
         % solving using usual x=A\b command
         dYt = JYt \-fYt;

         Yt_new  = Yt + alpha .* dYt;
         k   = k +1;
         error2 = norm (dYt);
         error3 = norm (MLSS_NonLin(Yt,Params,rho));
         error4 = norm (MLSS_NonLin(Yt_new,Params,rho));
         
         % Adaptive step size
         % this also can be used for MS4 task D
         if( error4 > error3  & alpha >1e-13)
             alpha = alpha / 10;
             if (alpha < 1e-13)
                 fprintf('ERROR : Newton is not converging');
                 break;
             end
         else
             Yt = Yt_new;
             error1 = error2;
             alpha = alpha * 1.01;
         end
    end

newton_sol = Yt;

end