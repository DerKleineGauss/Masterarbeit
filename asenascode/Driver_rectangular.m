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
Nx = 3; Ny = 3;
Nfaces = 4;

% Read in Mesh
Kx= 2;
Ky= 5;
[Nv, VX, VY, K, EToV, BCType] = rectangularGrid(Lq, Lr);

% Initialize solver and construct grid and metric
StartUp2D_rectangular;

% set up boundary conditions
BuildBCMaps2D;

B= functionB(x,y,a0, U,Lq,g,w,W0,n,delta, L_D);

% figure(4)
% plot3(x(:),y(:), B(:),'x')


[A] = PoissonIPDG2D(sparse(1:K*Np, 1:K*Np, B(:), K*Np, K*Np));

% set up Dirichlet boundary conditions
a= gamma*m*kB*Temp/(2*pi*hbar)^2;
b= epsilon*hbar/kB/Temp;
c= mu/kB/Temp;
% bed= Lr/2-abs(Fx(mapD)) < NODETOL;
% fh= @(k) a.*cos(Fy(mapD(bed)).*k).*log(1+exp(-b.*k.^2+c));
% f= integral(fh,-2*c,2*c,'ArrayValued', true);
bed= Lr/2-abs(Fx(mapD)) < NODETOL;
fh= @(k) a.*cos(Fy(mapD(bed)).*k).*log(1+exp(-b.*k.^2+c));
f= integral(fh,-2*c,2*c,'ArrayValued', true);

% figure(101)
% plot(Fy(mapD(bed)),f,'x')


rho_D= zeros(Nfp*Nfaces, K);

rho_D(mapD(bed))= f;
uD= rho_D;


% figure(2)
% bed2= abs(Lr/2+Fx(mapD)) < NODETOL;
% plot(Fy(mapD(bed2)), uD(mapD(bed2)),'x')


% set up Neumann boundary conditions
qN = zeros(Nfp*Nfaces, K);
% qN(mapN) = nx(mapN).*(pi*cos(pi*Fx(mapN)).*sin(pi*Fy(mapN))) + ...
%            ny(mapN).*(pi*sin(pi*Fx(mapN)).*cos(pi*Fy(mapN))) ;

% evaluate boundary condition contribution to rhs
Aqbc = PoissonIPDGbc2D(uD, qN);

% set up right hand side forcing
% rhs = -2*(pi^2)*sin(pi*x).*sin(pi*y);
rhs= 0;
rhs = -MassMatrix*(J.*rhs) + Aqbc;

% solve system
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
xlabel('$r / nm$')
ylabel('$q / nm$')
zlabel('$\rho(r,q)$')

figure(2)
rho_L= Z(abs(Y) <= NODETOL);
xl= X(abs(Y) <= NODETOL);
plot(xl,rho_L)
xlabel('$x / nm$')
ylabel('$n(x)$')

toc

tol=1e-14;
constant_r_vec = uniquetol(x,tol);
block_index = 0;
for constant_r = constant_r_vec
    block_index = block_index + 1;
    mask = find(abs(x-constant_r)<tol);
    N = length(mask);
end

