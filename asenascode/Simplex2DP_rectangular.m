function [P] = Simplex2DP_rectangular(r,s,i,j);

% function [P] = Simplex2DP(a,b,i,j);
% Purpose : Evaluate 2D orthonormal polynomial
%           on simplex at (a,b) of order (i,j).

h1 = JacobiP(r,0,0,i);
h2 = JacobiP(s,0,0,j);
P = h1.*h2;
return;
