function [f] = LvN_f(r,q,u)
% assume r,q,u have same dimensions

% alle Angaben nach Lukas paper
e = 1.60217662e-19;
m=9.1e-31;
hbar=6.626070040e-34/(2*pi);

L1_half = 3e-9; % nm
L2 = 2.5e-9; % nm
V_p = 0*r;
V_m = 0*r;
max_V = 0.1768*e;  % eV
V_p(abs(r+q/2)>=L1_half & abs(r+q/2)<=(L1_half+L2)) = max_V;
V_m(abs(r-q/2)>=L1_half & abs(r-q/2)<=(L1_half+L2)) = max_V;
% V = 0;
f = e*m/(hbar^2) * (V_p - V_m) .* u;

end
