function [B] = functionB(params, time)

r = params.x_interface;
q = params.y_interface;

a0 = params.constants.a0;
U = voltage(time, params.U, params.rampTime);
W0 = params.constants.W0;
n = params.constants.n;

Lq = params.Lq_scaled;
g = params.g_scaled;
w = params.w_scaled;
delta = params.delta_scaled;
L_D = params.L_D_scaled;

x = linspace(-params.Lr_scaled/2, params.Lr_scaled/2, 1000);
% plot(x,Potential(x,a0,U, g, w, L_D))

B= (Potential(r+0.5*q,a0,U, g, w, L_D)...
    -Potential(r-0.5*q,a0,U, g, w, L_D))...
    -1i*functionW(q, Lq/2-delta, W0, n, delta);

% B = 0*B;
end
