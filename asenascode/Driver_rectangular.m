setLatex;
close, clear()
clc

tic
% Driver script for solving the 2D Poisson equation
params = GlobalParams;
params.testing = true;

%% Umskalierung
epsilon= params.constants.e / params.constants.hbar;
gamma= sqrt(epsilon*2*params.constants.m/params.constants.hbar);
[params.Lr_scaled, params.Lq_scaled, params.L_D_scaled, params.w_scaled, params.g_scaled, params.delta_scaled] = params.scale(gamma);

% Polynomial order used for approximation, Grid size
params.Nx = 3;
params.Ny = 3;
params.Nfaces = 4;
params.Kx= 20;
params.Ky= 20;
% derived values
params.Npx = params.Nx + 1; params.Npy = params.Ny + 1; params.Nfp = params.Npx;
params.Np = params.Npx * params.Npy;
params.K = params.Ky * params.Kx;
% Read in Mesh
params = rectangularGrid(params);

% Initialize solver and construct grid and metric
params = StartUp2D_rectangular(params);

% Set up DFT matrix Phi
[Phi, params.p_DFT] = BuildPhi_y(params);
params.p_DFT = reshape(params.p_DFT, params.Np, params.K);
params.P_DFT = params.p_DFT(params.Fmask(:),:);    % face p's
% adapt boundarys
[params]=adaptBoundaryMap(params);
if(params.testing)
    figure(10)
    plot_sorted(params.x(:), params.p_DFT(:))
end

B= functionB(params);

if(params.testing)
    figure(4)
    plot3(params.x(:),params.y(:), B(:),'x')
end

[A] = PoissonIPDG2D_rectangular(params, sparse(1:params.K*params.Np, 1:params.K*params.Np, B(:), params.K*params.Np, params.K*params.Np));

%% boundary conditions
% set up Dirichlet boundary conditions
[a, b, c] = params.get_abc(gamma, epsilon);

% old bc's in space regime
% bed = Lr/2-abs(Fx(mapD)) < NODETOL;
% fh= @(k) a.*cos(Fy(mapD(bed))*k).*log(1+exp(-b.*k^2+c));
% f= integral(fh,-2*c,2*c,'ArrayValued', true);
%
% figure(101)
% plot_sorted(Fy(mapD(bed)),f)
%
% uD= zeros(Nfp*Nfaces, K);
% uD(mapD(bed))= f;

% new bc's in k regime
uD = rhs_dft(params, a,b,c, Phi);

% set up Neumann boundary conditions
qN = zeros(params.Nfp*params.Nfaces, params.K);

% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D_rectangular(params, uD, qN);

% set up right hand side forcing
% rhs = -2*(pi^2)*sin(pi*x).*sin(pi*y);
rhs= 0;
rhs = -params.MassMatrix*(params.J.*rhs) + Aqbc;

%% solve system
u_hat = (Phi*A*Phi')\rhs(:);
u = Phi' * u_hat;
u = reshape(u, params.Np, params.K);

% figure(1)
% %hack
% %xi = sqrt(hbar^2 / m / (a0*e));
% x = params.x; y = params.y;
% x = x/gamma*2*1e9;
% y = y/gamma*2*1e9;
% u = u*2;
% %
% x1= x(:);
% x2= y(:);
% m1= min(x1);
% m2= max(x1);
% dx= (m2-m1)/400;
% xlin= m1:dx:m2;
% m1= min(x2);
% m2= max(x2);
% dy= (m2-m1)/400;
% ylin= m1:dy:m2;
% [X, Y]= meshgrid(xlin, ylin);
% f= scatteredInterpolant(x(:), y(:),real(u(:)),'linear');
% %f= scatteredInterpolant(x(:), y(:),imag(u(:)),'linear');
% Z= f(X,Y);
% mesh(X,Y, f(X,Y))
% xlabel('$r / \xi$')
% ylabel('$q / \xi$')
% zlabel('$\rho(r,q)$')
% 
% figure(2)
% rho_L= Z(abs(Y) <= params.NODETOL);
% xl= X(abs(Y) <= params.NODETOL);
% plot(xl,rho_L)
% xlabel('$x / \xi$')
% ylabel('$n(x)$')

%eigenwerte = eigs(A)

toc
