function [P] = Simplex2DP(r,s,i,j);

% function [P] = Simplex2DP(r,s,i,j);
% Purpose : Evaluate 2D orthonormal polynomial
%           on cube at (r,s) of order (i,j).

h1 = JacobiP(r,0,0,i);
h2 = JacobiP(s,0,0,j);
P = h1.*h2;
return;
