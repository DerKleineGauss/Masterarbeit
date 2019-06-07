function plot_solution(params, u)

figure(1)
%hack
%xi = sqrt(hbar^2 / m / (a0*e));
x = params.x; y = params.y;
gamma = params.gamma;

x = x/gamma*2*1e9;
y = y/gamma*2*1e9;
% u = u*2;
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