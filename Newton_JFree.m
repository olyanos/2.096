function [newton_sol,k] = Newton_JFree(Yt,Gamma,Lambda,A,Params,rho, delta, eps)
% find the 

k = 0;
n = length(Yt);
alpha = 1;
fYt = MLSS_NonLin(Yt,Gamma,Lambda,A,Params,rho);
dYt = ones(n,1);
error1 = 1;
while ( norm (dYt) > 1e-3 & norm (fYt) > 1e-3)
         fYt = MLSS_NonLin(Yt,Gamma,Lambda,A,Params,rho); % function calculated @ (Yt);           
         % solving using Jacobian free GCR
         fhand = @(x) MLSS_NonLin(x,Gamma,Lambda,A,Params,rho);
         dYt = tgcr_(fhand, -fYt, delta, eps, 1000);

         Yt_new  = Yt + alpha .* dYt;
         k   = k +1;
         error2 = norm (dYt);
         error3 = norm (MLSS_NonLin(Yt,Gamma,Lambda,A,Params,rho));
         error4 = norm (MLSS_NonLin(Yt,Gamma,Lambda,A,Params,rho));
         
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
             alpha = alpha * 1.2;
         end
end

newton_sol = Yt;

end