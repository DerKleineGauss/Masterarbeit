function [W] = functionW(q, q0, W0, n, delta)

W= W0/delta^(2*n).*((q-q0).^(2*n).*(q >= q0 & q <= q0+delta)...
    + (q+q0).^(2*n).*(q >= -q0-delta & q <= -q0));

end

