function [y] = fermi_dirac_ft(k, q, mu,T, m, kb, hbar)
y = cos(k.*q).*fermi_dirac(k,mu,T, m, kb, hbar);
end