function [nx] = Normals1D(params)

% function [nx] = Normals1D
% Purpose : Compute outward pointing normals at elements faces

nx = zeros(params.Nfp*params.Nfaces, params.K); 

% Define outward normals
nx(1, :) = -1.0; nx(2, :) = 1.0;
return
