function [x] = Build_x_EXP(y,yw,P,Rx,Gamma,Lambda)
%BUILD_X_EXP     Build flow vector x

% Inputs:         1) y:      vector of nodal GDPs
%                 2) yw:     world GDP
%                 3) P:      nodal consumer price indices
%                 4) Rx:     bilateral trade resistances (exports)
%                 5) Gamma:  selection matrix
%                 6) Lambda: summation matrix
%
% Outputs:        1) x:      Flow vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = zeros(n*m,1);
% kk = 1;
% for ii = 1:n
%     for jj = 1:m
%         % We are dealing with the line connecting nodes ib and jb, where 
%         % ib = ii and jb = jj if ii>jj, else, jb = jj+1
%         ib = ii;
%         if ii <= jj
%             jb = jj + 1;
%         else
%             jb = jj;
%         end
%         x(kk) = (y(ib)*y(jb)/yw)*(P(ib)*P(jb))/Rx(kk);
%         
%         % Increment
%         kk = kk + 1;
%     end
% end

% Written via matrices
x = (1/yw)*(Lambda'*y).*(Gamma*y).*(Lambda'*P).*(Gamma*P)./Rx;
end

