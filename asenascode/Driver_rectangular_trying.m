setLatex;
close, clear()
clc

tic
% Driver script for solving the 2D Poisson equation
Globals2D;
testing = false;

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
Nx = 3; Ny = 3;
Npx = Nx + 1; Npy = Ny + 1;
Nfaces = 4;

% Read in Mesh
Kx= 50;
Ky= 50;
[Nv, VX, VY, K, EToV, BCType] = rectangularGrid(Lq, Lr);

% Initialize solver and construct grid and metric
StartUp2D_rectangular;

% Set up DFT matrix Phi
[Phi, p_DFT] = BuildPhi_y();
p_DFT = reshape(p_DFT, Np, K);
P_DFT=p_DFT(Fmask(:),:);    % face p's
% adapt boundarys
% adaptBoundaryMap(Lr)
if(testing)
    figure(10)
    plot_sorted(x(:), p_DFT(:))
end

B= functionB(x,y,a0, U,Lq,g,w,W0,n,delta, L_D);

if(testing)
    figure(4)
    plot3(x(:),y(:), B(:),'x')
end


[A] = PoissonIPDG2D_rectangular(sparse(1:K*Np, 1:K*Np, B(:), K*Np, K*Np));

%% boundary conditions
% set up Dirichlet boundary conditions
a= gamma*m*kB*Temp/(2*pi*hbar)^2;
b= epsilon*hbar/kB/Temp;
c= mu/kB/Temp;

% old bc's in space regime
bed = Lr/2-abs(Fx(mapD)) < NODETOL;
fh= @(k) a.*cos(Fy(mapD(bed)).*k).*log(1+exp(-b.*k.^2+c));
f= integral(fh,-2*c,2*c,'ArrayValued', true);

if (testing)
    figure(101);
    plot_sorted(Fy(mapD(bed)),f);
end

uD= zeros(Nfp*Nfaces, K);
uD(mapD(bed))= f;

% % new bc's in k regime
% uD = rhs_dft(a,b,c, Lr, Phi);

% set up Neumann boundary conditions
qN = zeros(Nfp*Nfaces, K);

% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D_rectangular(uD, qN);
% Aqbc = Phi * Aqbc(:);
% Aqbc = reshape(Aqbc, Np, K);

% set up right hand side forcing
% rhs = -2*(pi^2)*sin(pi*x).*sin(pi*y);
rhs= 0;
rhs = -MassMatrix*(J.*rhs) + Aqbc;

%% solve system
% u_hat = (Phi*A*Phi')\rhs(:);
% u = Phi' * u_hat;

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

%eigenwerte = eigs(A)

toc
