setLatex;
close, clear()
clc

tic
% Driver script for solving the 2D Poisson equation
Globals2D;


%% Konstanten
h= 5e-9;
Lr= 106e-9/2;
Lq= 152e-9/2;
L_D= 30e-9/2;
w= 5e-9;
g= 5.5e-9;

hbar= 1.054571800e-34;
e= 1.6021766208e-19;
U= 0;
delta= 0.2*Lq;
W0= 1;
n= 1;
a0= 0.1768;

m0= 9.10938356e-31;
m= 0.063*m0;

kB= 1.38064852e-23;
Temp= 300;
mu= 0.0467*e;

kappa= 1e8;
%% Umskalierung
epsilon= e/hbar;
gamma= sqrt(epsilon*2*m/hbar);

h= h*gamma;
epsh= e-12;
Lr= floor(Lr*gamma);
Lq= floor(Lq*gamma);
L_D= L_D*gamma;
w= w*gamma;
g= g*gamma;
delta= delta*gamma;



% Polynomial order used for approximation 
Nx = 1; Ny = 1;
Npx = Nx + 1; Npy = Ny + 1;
Nfaces = 4;

% Read in Mesh
Kx= 2;
Ky= 4;
[Nv, VX, VY, K, EToV, BCType] = rectangularGrid(Lq, Lr);

% Initialize solver and construct grid and metric
StartUp2D_rectangular;

% Set up DFT matrix Phi
[Phi, p_DFT] = BuildPhi_y(x, y, Kx, Ky, Npx, Npy);
p_DFT = reshape(p_DFT, Np, K);
P_DFT=p_DFT(Fmask(:),:);    % face p's
% figure(10)
% plot_sorted(x(:), p_DFT)

% set up boundary conditions
BuildBCMaps2D;

B= functionB(x,y,a0, U,Lq,g,w,W0,n,delta, L_D);

% figure(4)
% plot3(x(:),y(:), B(:),'x')


[A] = PoissonIPDG2D_rectangular(sparse(1:K*Np, 1:K*Np, B(:), K*Np, K*Np));

%% boundary conditions
% set up Dirichlet boundary conditions
a= gamma*m*kB*Temp/(2*pi*hbar)^2;
b= epsilon*hbar/kB/Temp;
c= mu/kB/Temp;

map_left_v = Lr/2 + x < NODETOL;
map_right_v = Lr/2 - x < NODETOL;
map_left_f = Lr/2 + Fx < NODETOL;
map_right_f = Lr/2 - Fx < NODETOL;

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
right_pos = find(Lr/2-VX < NODETOL & VY <= 0);
left_pos  = find(Lr/2+VX < NODETOL & VY <= 0);
right_neg = find(Lr/2-VX < NODETOL & VY >= 0);
left_neg  = find(Lr/2+VX < NODETOL & VY >= 0);
plot(VX,VY,'k.',VX(right_pos), VY(right_pos), 'ro', VX(left_pos), VY(left_pos), 'bo', ...
    VX(right_neg), VY(right_neg), 'yx', VX(left_neg), VY(left_neg), 'gx');
CorrectBCTable(left_neg, Neuman);
CorrectBCTable(right_pos, Neuman);
BuildBCMaps2D;
plot(Fx(:),P_DFT(:),'k.',Fx(mapN),P_DFT(mapN),'ro', Fx(mapD),P_DFT(mapD),'gx');
legend('all face points','Neumann face points', 'Dirichlet face points', 'Location','north')

fh_left = @(k) a.*cos(y(map_left_v )*k).*log(1+exp(-b.*k^2+c));
fh_right= @(k) a.*cos(y(map_right_v)*k).*log(1+exp(-b.*k^2+c));
% uD_v caontains values in "volume speak" generated with x and y instead Fx
% and Fy (as is the case for the later defined uD_complete)
uD_v = zeros(Np, K);
uD_v(map_left_v) = integral(fh_left,-2*c,2*c,'ArrayValued', true);
uD_v(map_right_v) = integral(fh_right,-2*c,2*c,'ArrayValued', true);
% left and right sides must have same y values -> check here
assert ( norm(uD_v(map_left_v) - uD_v(map_right_v)) < NODETOL)
figure(111);
plot_sorted(y(map_left_v), uD_v(map_left_v));
% confirm we do not lose map control
assert(norm(find(Phi*uD_v(:))-find(uD_v)) < NODETOL)
% apply DFT to rhs
uD_v = reshape(Phi*uD_v(:), Np, K);
% go to needed definitions for each face as also done in StartUp2D for Fx
uD_complete = uD_v(Fmask(:), :);
uD = zeros(Nfp*Nfaces, K);

% set up Neumann boundary conditions
qN = zeros(Nfp*Nfaces, K);

% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D_rectangular(uD, qN);

% set up right hand side forcing
% rhs = -2*(pi^2)*sin(pi*x).*sin(pi*y);
rhs= 0;
rhs = -MassMatrix*(J.*rhs) + Aqbc;

%% solve system
u = A\rhs(:);
u = reshape(u, Np, K);

figure(1)
%hack
%xi = sqrt(hbar^2 / m / (a0*e));
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
rho_L= Z(abs(Y) <= NODETOL);
xl= X(abs(Y) <= NODETOL);
plot(xl,rho_L)
xlabel('$x / \xi$')
ylabel('$n(x)$')

eigenwerte = eigs(A)

toc

