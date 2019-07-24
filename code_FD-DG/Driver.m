clear;
close all;
tic

% Driver script for solving the LVN equation

%% set parameters
params = getDefaultParams();


%% Generate simple mesh
xmin = -params.Lr_scaled/2;
xmax = +params.Lr_scaled/2;
[params.VX, params.EToV] = MeshGen1D(xmin, xmax, params.K);
params.hx = (params.VX(1,2)-params.VX(1,1)) ;

ymin = -params.Lq_scaled/2;
ymax = +params.Lq_scaled/2;
params.y_interface = linspace(ymin, ymax, params.Npy);
params.hy = params.y_interface(2)-params.y_interface(1);
params.y = linspace(ymin+params.hy/2, ymax-params.hy/2, params.Ny);

%% Initialize solver and construct grid and metric
params = StartUp1D(params);

% adapt grid to 2D implementation
[params.y_interface, params.x_interface] = meshgrid(params.y_interface,params.x);
[params.y, params.x] = meshgrid(params.y, params.x);

%% get hopefully final equilibrium solution in v
B= functionB(params, params.FinalTime);
[A, rhs, R] = LVN_systemmatrix(params, B);
% Solve Problem
v_final = A\rhs;
v_final = reshape(v_final,params.Np*params.K, params.Ny);

%% initial and timestepping
if (params.timestepping)
    % get initial equilibrium solution in v
    % set potential
    B= functionB(params, -1);
    [A, rhs, R] = LVN_systemmatrix(params, B);
    % Solve Problem
    v = A\rhs;
    v = reshape(v,params.Np*params.K, params.Ny);

    % solve time stepping and transform v -> u
    v = timeStepping(params, v(:), v_final);
    v = reshape(v,params.Np*params.K, params.Ny);
    u = strangeMatrixMultiplication(R , v);
    u = reshape(u,params.Np*params.K, params.Ny);
    plot_solution(params, u, true);
    plot_solution(params, u, false);
else
    u = strangeMatrixMultiplication(R , v_final);
    u = reshape(u,params.Np*params.K, params.Ny);
    plot_solution(params, u, true);
    plot_solution(params, u, false);
    fid = figure('name', 'Dichte');
    plotDensity_FVDG(v_final, params, R, 300, fid, params.FinalTime, 'b');
end

toc
