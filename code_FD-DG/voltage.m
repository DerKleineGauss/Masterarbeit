function [U] = voltage(t, U0, rampTime)
if (t<0)
    U = 0;
elseif (t<rampTime)
    U = t/rampTime * U0;
else
    U = U0;
end
end