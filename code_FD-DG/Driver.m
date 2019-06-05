clear;
close all;
tic

% Driver script for solving the LVN equation

params = GlobalParams;
params.testing = true;

% Rescaling
epsilon= params.constants.e / params.constants.hbar;
gamma= sqrt(epsilon*2*params.constants.m/params.constants.hbar);
[params.Lr_scaled, params.Lq_scaled, params.L_D_scaled, params.w_scaled, params.g_scaled, params.delta_scaled] = params.scale(gamma);

% Order of polymomials used for approximation 
params.N = 2;
params.K = 6;
params.Npy = 7;

% Generate simple mesh
xmin = -params.Lr_scaled/2;
xmax = +params.Lr_scaled/2;
[params.VX, params.EToV] = MeshGen1D(xmin, xmax, params.K);

ymin = -params.Lq_scaled/2;
ymax = +params.Lq_scaled/2;
params.y = linspace(ymin, ymax, params.Npy);
params.hy = params.y(2)-params.y(1);

% Initialize solver and construct grid and metric
params = StartUp1D(params);

% adapt grid to 2D implementation
[params.y, params.x] = meshgrid(params.y,params.x);

% set potential
B= functionB(params);
if(params.testing)
    figure(4)
    plot3(params.x(:),params.y(:), real(B(:)),'x')
    figure(5)
    plot3(params.x(:),params.y(:), imag(B(:)),'x')
end

[A] = LVN_systemmatrix(params, B);

% Set initial conditions
% u = sin(x);

% Solve Problem



