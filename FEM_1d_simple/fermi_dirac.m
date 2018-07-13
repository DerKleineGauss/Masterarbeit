function [f] = fermi_dirac(k,mu,T, m, kb, hbar)
f = m*kb*T/(2*pi*hbar^2) * log(1+exp(-1/(kb*T)*(hbar^2*k.^2/(2*m) - mu)));
%f = log(1+exp(-1/(kb*T)*(hbar^2*k.^2/(2*m) - mu)));
%f = 1./(1+exp(1/(kb*T)*(hbar^2*k.^2/(2*m) - mu)));
end