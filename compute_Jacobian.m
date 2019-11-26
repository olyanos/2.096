function [JYt] = compute_Jacobian(Yt,Gamma,Lambda,A,Params,rho,h);

n = length(Yt);
JYt = zeros(n,n);
fYt = MLSS_NonLin(Yt,Params,rho);
    for k = 1:n
        Yt_new = Yt;
        Yt_new(k) = Yt_new(k)+h;
        fYt_new = MLSS_NonLin(Yt_new,Params,rho);
        JYt(:,k) = (fYt_new-fYt)./h;
    end
    
end