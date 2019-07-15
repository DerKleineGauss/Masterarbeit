function [B] = functionB(params, time)

r = params.x_interface;
q = params.y_interface;

% x = linspace(-params.Lr_scaled/2, params.Lr_scaled/2, 1000);
% plot(x,Potential(x,a0,U, g, w, L_D))

if (params.withCAP)
    B= (Potential(r+0.5*q, time, params)...
       -Potential(r-0.5*q, time, params))...
       -1i*functionW(q, params);
else
    B= (Potential(r+0.5*q, time, params)...
       -Potential(r-0.5*q, time, params));
end

% B = 0*B;
end
