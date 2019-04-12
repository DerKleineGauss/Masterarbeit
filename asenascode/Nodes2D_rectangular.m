function [r,s] = Nodes2D_rectangular(Nx, Ny);

% function [x,y] = Nodes2D_rectangular(Nx, Ny);
% Purpose  : Compute (x,y) nodes in reference cube for polynomial of order Nx,Ny

x_linear = JacobiGL(0,0,Nx);
y_linear = JacobiGL(0,0,Ny);

[r,s] = meshgrid(x_linear, y_linear);
s = s(:);
r = r(:);

return;
