function [B] = B (r,q)
B = potential(r + q/2) - conj(potential(r - q/2));
end