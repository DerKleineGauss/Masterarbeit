function [Dr,Ds] = Dmatrices2D_rectangular(params)

% function [Dr,Ds] = Dmatrices2D(N,r,s,V)
% Purpose : Initialize the (r,s) differentiation matrices
%	    on the simplex, evaluated at (r,s) at order N

Npx = params.Npx;
Npy = params.Npy;
r = params.r;
s = params.s;
V = params.V;

[Vr, Vs] = GradVandermonde2D_rectangular(params);
Dr = Vr/V; Ds = Vs/V;
return;
