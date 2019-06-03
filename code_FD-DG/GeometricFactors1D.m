function [rx,J] = GeometricFactors1D(params)

% function [rx,J] = GeometricFactors1D(x,Dr)
% Purpose  : Compute the metric elements for the local mappings of the 1D elements     

xr  = params.Dr*params.x; J = xr; rx = 1./J; 
return
