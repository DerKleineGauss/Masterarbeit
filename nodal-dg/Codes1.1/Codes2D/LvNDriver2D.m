% Driver script for solving the 2D Poisson equation
Globals2D;
GlobalsLvN;

% Polynomial order used for approximation
N = 5;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('block2.neu');

% Initialize solver and construct grid and metric
StartUp2D;

% set up boundary conditions
BuildBCMaps2D;

% set up right hand side for homogeneous Poisson
[A,M] = Poisson2D(); % Setup using PoissonRHS2D.m
[A,M] = PoissonIPDG2D(); % Setup using PoissonIPDG2D.m

% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);
mapSides =  find(abs(abs(Fx)-max(Fx(:))) < 1e-12);
mapTopBot = find(abs(abs(Fy)-max(Fy(:))) < 1e-12 & abs(abs(Fx)-max(Fx(:))) > 1e-12);
q = Fy(mapSides) * L_q; % 1d real coordinates at the two sides left and right

% Calculate mu
mu = newtonRaphson(@nullstellenSucheMu, 1.5*e);

f_hut = zeros(size(q,1),1);
f_hut_trapez = zeros(size(q,1),1);
fermDiracFt = @(k_value, q_value) fermi_dirac_ft(k_value, q_value, mu);
upper_k = sqrt(2*m*mu/hbar/hbar);
safety_k = 1.3;
nr_steps_k = 5e6;
steplength_k = safety_k*upper_k/nr_steps_k;

k = linspace(0 , safety_k*upper_k , nr_steps_k);
for i=1:size(q,1)
%     f_hut(i) = 2/(2*pi)*integral( @(k)fermDiracFt(k, q(i)),0,2e10);
%     f_hut(N_q-i+1) = f_hut(i);
    y_temp = fermDiracFt(k, q(i));
    f_hut_trapez(i) = steplength_k * 2/(2*pi)*trapz(y_temp);
end
clear k y_temp;

uD(mapSides) = f_hut_trapez;
uD(mapTopBot) = 0;

% evaluate boundary condition contribution to rhs
qN = zeros(Nfp*Nfaces, K);
Aqbc = PoissonIPDGbc2D(uD, qN);

% set up right hand side forcing
rhs = -2*(pi^2)*sin(pi*x).*sin(pi*y);
rhs = -MassMatrix*(J.*rhs) + Aqbc;

% solve system
u = A\rhs(:);
u = reshape(u, Np, K);
