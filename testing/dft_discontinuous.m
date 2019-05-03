Np = 2;
K = 500;
a = 3;

%% discontinuous case
x = zeros(K*Np,1);
u = zeros(K*Np,1);
xmax = 1 *2*pi;
deltax = xmax / K / (Np-1);
for k=1:K
    xL = (k-1)*xmax/K;
    xR = k*xmax/K;
    x((k-1)*Np + 1 : k*Np) = xL: deltax :xR;
    u((k-1)*Np + 1 : k*Np) = sin(a*(xL+xR)/2);
    %u((k-1)*Np+1) = sin(a*xL);
    %temp = a*cos(a*xL);
    %u((k-1)*Np + 2 : k*Np) = sin(a*xL) + temp*( x((k-1)*Np + 2 : k*Np) - xL);
end
%u=sin(a*x);
% figure(4)
% plot(x,u);
% figure(5)
% plot(x,sin(a*x));
y1 = fft(u);
plotdft(y1,xmax,1)

% try less k points
N_k = K;%Np*K-K;  % number of points in k space
k = 2*pi/xmax * (0:N_k-1);
[X_mesh, K_mesh] = meshgrid(x,k);
Phi = exp(-1i*K_mesh.*X_mesh) / sqrt(length(u));
y3 = Phi*u;
f = k;
plotdft_with_f(y3, f, 2);
figure(8)
subplot(2,1,1)
plot(x, (Phi'*y3));
title('\Phi^T \Phi u')

subplot(2,1,2)
plot(x,u);
title('u')


% try factor 1/2 where x_k=x_k+1 and same k value
% N_k = length(u);
% kmax = 2*pi/xmax;
% k=zeros(N_k,1);
% deltak = kmax / K / (Np-1);
% for e=1:K
%     kL = kmax*(e-1)/K;
%     kR = kmax*e/K;
%     k((e-1)*Np + 1 : e*Np) = kL: deltak :kR;
% end

%% continuous case
xmax=2*pi*10;
%N = Np*K-K;
N=1200;
x2 = 0:xmax/N:xmax-(xmax/N);
u2=sin(a*x2);
y2 = fft(u2);
plotdft(y2,xmax,5)


function plotdft(y, xmax, fignr)
    figure(fignr);
    m = abs(y);                               % Magnitude
    y(m<1e-6) = 0;
    p = unwrap(angle(y));                     % Phase

    f = (0:length(y)-1)/xmax*2*pi;        % Frequency vector

    subplot(2,1,1)
    plot(f,m)
    title('Magnitude')

    subplot(2,1,2)
    plot(f,p*180/pi)
    title('Phase')
end

function plotdft_with_f(y, f, fignr)
    figure(fignr);
    m = abs(y);                               % Magnitude
    y(m<1e-6) = 0;
    p = unwrap(angle(y));                     % Phase

    subplot(2,1,1)
    plot(f,m)
    title('Magnitude')

    subplot(2,1,2)
    plot(f,p*180/pi)
    title('Phase')
end