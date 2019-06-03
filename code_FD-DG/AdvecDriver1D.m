% Driver script for solving the 1D advection equations
params = GlobalParams;

% Order of polymomials used for approximation 
params.N = 2;
params.K = 10;

% Generate simple mesh
[params.VX, params.EToV] = MeshGen1D(0.0,2.0,params.K);

% Initialize solver and construct grid and metric
params = StartUp1D(params);

% Set initial conditions
u = sin(params.x);

% Solve Problem
FinalTime = 40;
figure(1);

[u] = Advec1D(u,FinalTime, params);


