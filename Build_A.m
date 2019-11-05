function [A] = Build_A(n)
% BUILD_A    Build matrix A which converts from x to q via q = Ax

% Inputs:    1) n: number of nodes
% Outputs:   1) A: conversion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
m = n-1;
A = zeros(n*m,n*m);

% Define original line vector
bv = zeros(n,m);
for gg = 1:n
    bv(gg,:) = (1:m) + (gg-1)*m;
end

% loop over nodes and lines
kk = 1;
for ii = 1:n
    for jj = 1:m
        if ii <= jj
            jb = jj + 1;
            if ii <= jj
                A(kk,bv(jb,ii)) = 1;
            else
                A(kk,bv(jb,ii-1)) = 1;
            end
        else
            jb = jj;
            if ii <= jj
                A(kk,bv(jb,ii)) = 1;
            else
                A(kk,bv(jb,ii-1)) = 1;
            end
        end
        % Increment
        kk = kk + 1;
    end
end
end

