function [Dr] = Dmatrix1D(params)

% function [Dr] = Dmatrix1D(N,r,V)
% Purpose : Initialize the (r) differentiation matrices on the interval,
%	        evaluated at (r) at order N

Vr = GradVandermonde1D(params.N, params.r);
Dr = Vr/params.V;
return
