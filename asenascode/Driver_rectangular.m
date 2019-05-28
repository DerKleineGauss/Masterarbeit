setLatex;
close, clear()
clc
opengl('save', 'software')

tic
% Driver script for solving the 2D Poisson equation
params = GlobalParams;
params.testing = false;

% Rescaling
epsilon= params.constants.e / params.constants.hbar;
gamma= sqrt(epsilon*2*params.constants.m/params.constants.hbar);
[params.Lr_scaled, params.Lq_scaled, params.L_D_scaled, params.w_scaled, params.g_scaled, params.delta_scaled] = params.scale(gamma);

% Polynomial order used for approximation, Grid size
params.Nx = 2;
params.Ny = 2;
params.Nfaces = 4;
params.Kx= 2;
params.Ky= 2;
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
[params] = BuildBCMaps2D(params);

if(params.testing)
    figure(10)
    plot(params.x(:), params.p_DFT(:), 'x')
end

B= functionB(params);

if(params.testing)
    figure(4)
    plot3(params.x(:),params.y(:), B(:),'x')
end

[A, D] = PoissonIPDG2D_rectangular(params, sparse(1:params.K*params.Np, 1:params.K*params.Np, B(:), params.K*params.Np, params.K*params.Np));

[D_l_pos, D_l_neg, D_r_pos, D_r_neg] = SplitDmatrix(params, D, Phi);
%% boundary conditions
% set up Dirichlet boundary conditions
[a, b, c] = params.get_abc(gamma, epsilon);

% % old bc's in space regime
% bed = params.Lr_scaled/2-abs(params.x(params.vmapD)) < params.NODETOL;
% fh= @(k) a.*cos(params.y(params.vmapD(bed))*k).*log(1+exp(-b.*k^2+c));
% f= integral(fh,-2*c,2*c,'ArrayValued', true);
% 
% if (params.testing)
%     figure(101)
%     plot_sorted(params.y(params.vmapD(bed)),f)
% end
% 
% uD= zeros(params.Np, params.K);
% uD(params.vmapD(bed))= f;

% new bc's in k regime
rhs = rhs_dft(params, a,b,c, Phi, D_l_pos, D_r_neg, D_l_neg, D_r_pos);

% set up Neumann boundary conditions
qN = zeros(params.Nfp*params.Nfaces, params.K);

% adapt left hand side to match fourier transform and subtract unknown
% dirichlet boundary conditions
% A = A*Phi' - D_l_neg - D_r_pos;
A = A*Phi';

%% solve system
u_hat = A\rhs;
% Phi*u =: u_hat  -->  u=Phi' * u_hat
u = Phi' * u_hat;
u = reshape(u, params.Np, params.K);

figure(1)
%hack
%xi = sqrt(hbar^2 / m / (a0*e));
x = params.x; y = params.y;
x = x/gamma*2*1e9;
y = y/gamma*2*1e9;
u = u*2;
%
x1= x(:);
x2= y(:);
m1= min(x1);
m2= max(x1);
dx= (m2-m1)/400;
xlin= m1:dx:m2;
m1= min(x2);
m2= max(x2);
dy= (m2-m1)/400;
ylin= m1:dy:m2;
[X, Y]= meshgrid(xlin, ylin);
f= scatteredInterpolant(x(:), y(:),real(u(:)),'linear');
%f= scatteredInterpolant(x(:), y(:),imag(u(:)),'linear');
Z= f(X,Y);
mesh(X,Y, f(X,Y))
xlabel('$r / \xi$')
ylabel('$q / \xi$')
zlabel('$\rho(r,q)$')

figure(2)
rho_L= Z(abs(Y) <= params.NODETOL);
xl= X(abs(Y) <= params.NODETOL);
plot(xl,rho_L)
xlabel('$x / \xi$')
ylabel('$n(x)$')

%eigenwerte = eigs(A)

toc
