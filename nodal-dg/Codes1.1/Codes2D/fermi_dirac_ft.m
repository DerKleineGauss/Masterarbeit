function [y] = fermi_dirac_ft(k, q, mu)
y = cos(k.*q).*fermi_dirac(k,mu);
end