function [B] = B (r,q, e)
B = potential(r + q/2, e) - conj(potential(r - q/2, e));
end