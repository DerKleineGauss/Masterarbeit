function [W] = functionW(q, params)

W0 = params.constants.W0;
n = params.constants.n;
Lq = params.Lq_scaled;
delta = params.delta_scaled;
q0 = Lq/2-delta;

W= W0/delta^(2*n).*((q-q0).^(2*n).*(q >= q0 & q <= q0+delta)...
    + (q+q0).^(2*n).*(q >= -q0-delta & q <= -q0));

end

