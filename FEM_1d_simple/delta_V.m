function [dV] = delta_V (r,q,delta_q)
dV = potential(r + q + 0.5*delta_q) - conj(potential(r - 0.5*q - 0.5*delta_q));
end