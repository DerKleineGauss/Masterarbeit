function [rx,sx,ry,sy,J] = GeometricFactors2D(params)

% function [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds)
% Purpose  : Compute the metric elements for the local mappings of the elements

x = params.x;
y = params.y;
Dr = params.Dr;
Ds = params.Ds;
% Calculate geometric factors
xr = Dr*x; xs = Ds*x; yr = Dr*y; ys = Ds*y; J = -xs.*yr + xr.*ys;
rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
return;
