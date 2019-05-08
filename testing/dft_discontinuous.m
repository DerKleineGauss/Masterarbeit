clear(); close all;

Np = 2;
K = 50;
a = 10;

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
plotdft(y1,xmax,'discontinuous #k = #x')

%% try less k points
N_k = Np*K-K;  % number of points in k space
% N_k = K;  % number of points in k space
k = 2*pi/xmax * (0:N_k-1);
[X_mesh, K_mesh] = meshgrid(x,k);
Phi = exp(-1i*K_mesh.*X_mesh) / sqrt(length(u));
y3 = Phi*u;
f = k;
plotdft_with_f(y3, f, 'discontinuous #k < #x');

figure(8)
subplot(2,1,1)
plot(x, (Phi'*y3));
title('\Phi^T \Phi u')
subplot(2,1,2)
plot(x,u);
title('u')

hopefully_id = Phi'*Phi;
hopefully_id(abs(hopefully_id)<1e-12)=0;
display(['||Phi^T Phi - id|| = ', num2str(normest(hopefully_id - eye(size(hopefully_id))))]);


%% continuous case
N = Np*K;
% N=1200;
x2 = 0:xmax/N:xmax-(xmax/N);
u2=sin(a*x2);
y2 = fft(u2);
plotdft(y2,xmax,'continuous case')


function plotdft(y, xmax, fignr)
    y = fftshift(y);
    n = length(y);
    figure('name', fignr, 'NumberTitle','off');
    m = abs(y);                               % Magnitude
    y(m<1e-6) = 0;
    p = unwrap(angle(y));                     % Phase

    f = (-n/2:n/2-1)/xmax*2*pi;        % Frequency vector
    %f = fftshift(f);

    subplot(2,1,1)
    plot(f,m)
    title('Magnitude')

    subplot(2,1,2)
    plot(f,p*180/pi)
    title('Phase')
end

function plotdft_with_f(y, f, fignr)
    f = fftshift(f)-f(length(f)/2+1);
    figure('name', fignr, 'NumberTitle','off');
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