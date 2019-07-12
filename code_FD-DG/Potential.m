function [V] = Potential(x, time, params)

g = params.g_scaled;
a = params.constants.a0;
U = voltage(time, params.U, params.rampTime);
w = params.w_scaled;
L_D = params.L_D_scaled;

V= a.*(heaviside(x+w/2+g)-heaviside(x-w/2+g)+heaviside(x+w/2-g)-heaviside(x-w/2-g))...
    -U.*heaviside(x-L_D/2).*(-x./L_D+0.5)-U.*heaviside(x+L_D/2).*(x./L_D+0.5);

end

