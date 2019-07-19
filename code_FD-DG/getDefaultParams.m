function [params] = getDefaultParams()

params = GlobalParams;

% result and plot parameters
params.movieName = 'results/Tfinal_15000_N3_K220.avi';
params.makeMovie = true;
params.plot_logarithmic = true;
params.saveResults = false;

% Order of polymomials used for approximation (x direction)
params.N = 1;
% Number of cells for DG discretization (x direction)
params.K = 100;
% Number of cells (y direction)
params.Ny = 100;
% Number of interfaces (y direction)
params.Npy = params.Ny+  1;
% Voltage
params.U = +0.0;
params.rampTime = 0.2;
% final time
params.FinalTime = 150;
% time stepping mode ( one of rk4 , rk2ssp , )
params.mode = 'rk2ssp';
% do timestepping?
params.timestepping = true;

if (params.K * params.Ny * (params.N+1) > 150)
    params.testing = false;
else
    params.testing = true;
end

% BC related params
params.withCAP = false;
params.constants.mu_l = newtonRaphson(@nullstellenSucheMu, 1.5*params.constants.e, params);
% params.constants.E_c = -params.U*params.constants.e;
params.constants.mu_r = newtonRaphson(@nullstellenSucheMu, 1.5*params.constants.e, params);

% Rescaling
params.epsilon= params.constants.e / params.constants.hbar;
params.gamma= sqrt(params.epsilon*2*params.constants.m/params.constants.hbar);
[params.Lr_scaled, params.Lq_scaled, params.L_D_scaled, params.w_scaled, params.g_scaled, params.delta_scaled] = params.scale(params.gamma);
