clear;
close all;
tic

% Driver script for solving the LVN equation

params = GlobalParams;
params.testing = true;

% Rescaling
params.epsilon= params.constants.e / params.constants.hbar;
params.gamma= sqrt(params.epsilon*2*params.constants.m/params.constants.hbar);
[params.Lr_scaled, params.Lq_scaled, params.L_D_scaled, params.w_scaled, params.g_scaled, params.delta_scaled] = params.scale(params.gamma);

% Order of polymomials used for approximation (x direction)
params.N = 2;
% Number of cells for DG discretization (x direction)
params.K = 40;
% Number of cells (y direction)
params.Ny = 50;
% Number of interfaces (y direction)
params.Npy = params.Ny+  1;
% Voltage
% params.U = 0.25;
params.U = 0.25;
params.rampTime = 2;
% final time
params.FinalTime = 20;

if (params.K * params.Ny * params.N > 150)
    params.testing = false;
else
    params.testing = true;
end

% Generate simple mesh
xmin = -params.Lr_scaled/2;
xmax = +params.Lr_scaled/2;
[params.VX, params.EToV] = MeshGen1D(xmin, xmax, params.K);
params.hx = (params.VX(1,2)-params.VX(1,1)) / params.N;

ymin = -params.Lq_scaled/2;
ymax = +params.Lq_scaled/2;
params.y_interface = linspace(ymin, ymax, params.Npy);
params.hy = params.y_interface(2)-params.y_interface(1);
params.y = linspace(ymin+params.hy/2, ymax-params.hy/2, params.Ny);

% Initialize solver and construct grid and metric
params = StartUp1D(params);

% adapt grid to 2D implementation
[params.y_interface, params.x_interface] = meshgrid(params.y_interface,params.x);
[params.y, params.x] = meshgrid(params.y, params.x);

%% get initial equilibrium solution in v
% set potential
B= functionB(params, -1);
[A, rhs, R] = LVN_systemmatrix(params, B);
% Solve Problem
v = A\rhs;
v = reshape(v,params.Np*params.K, params.Ny);

%% solve time stepping and transform v -> u
v = timeStepping(params, v(:));
v = reshape(v,params.Np*params.K, params.Ny);
u = strangeMatrixMultiplication(R , v);
u = reshape(u,params.Np*params.K, params.Ny);
plot_solution(params, u);
% plot(params.x, u);

toc
