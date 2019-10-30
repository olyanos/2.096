function [Rx] = Build_BLT_Rx(gam1,gam2,d,nu,Y,rho,n,m)
%BUILD_BLT_RX    % Build the Bilateral Trade Resistance vector (exports)

% Inputs:       1) gam1: sensitivity vector (trade value) - per line!
%               2) gam2: sensitivity vector (currency) - per line!
%               3) d:    distance vector - per line!
%               4) nu:   commodity value - per line!
%               5) Y:    vector of nodal currencies - per node!
%               6) rho:  parameter scaling currencies
%               7) n:    number of nodes
%               8) m:    n-1
% Outputs:      1) Rx:   vector of BLT for exports - per line!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rx = zeros(n*m,1);
kk = 1;
for ii = 1:n      % loop over each node
    for jj = 1:m  % loop over each line
        % We are dealing with the line connecting nodes ib and jb, where 
        % ib = ii and jb = jj if ii>jj, else, jb = jj+1
        ib = ii;
        if ii <= jj
            jb = jj + 1;
        else
            jb = jj;
        end
        Rx(kk) = (gam1(kk)/nu(kk))*d(kk)*exp(rho*gam2(kk)*(Y(ib) - Y(jb)));
        
        % Increment
        kk = kk+1;
    end
end
end

