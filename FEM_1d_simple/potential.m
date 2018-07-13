function [V] = potential(x,e)
% alle Angaben nach Lukas paper
V = 0;
L1_half = 3e-9; % nm
L2 = 2.5e-9; % nm
max_V = 0.1768*e;  % eV
if abs(x)>=L1_half && abs(x)<=(L1_half+L2)
    V = max_V;
end
end