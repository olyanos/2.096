function [E,Lambda,Gamma] = Build_ELG(n)
%BUILD_ELG      This function builds the special incidence matrix
%               E associated with the system, along with the 

% Inputs:       1) n:      number of nodes
% Ouputs:       1) E:      special incidence matrix
%               2) Lambda: summation matrix
%               3) Gamma:  selection matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each of the n nodes has m=n-1 connections to the other nodes
m = n-1;

E = zeros(m*n,n);
for ii = 1:n
    em = -eye(m);
    ov = ones(m,1);
    am = [em(:,1:ii-1) ov em(:,ii:end)];
    E(((ii-1)*m)+(1:m),:) = am;
end

% Construct Lambda and Gamma
Lambda = 0.5*(E + abs(E))';
Gamma  = 0.5*(abs(E) - E);
end

