% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation 
N = 2;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0,2.0,10);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
u = sin(x);

% Solve Problem
FinalTime = 40;
figure(1);

[u] = Advec1D(u,FinalTime);


