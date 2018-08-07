function mu = getFermiLevel(n_dop, m_eff, T)

q   = 1.602E-19;
h   = 6.626E-34;
hb  = h/2/pi;
kB  = 1.38E-23;
m0  = 9.109E-31;

E_start = 0;
E_end   = 0.8;
dE      = 1E-4;
E       = E_start : dE : E_end;
mu      = (min(E)-0.2):dE:0.8;
dos_e   = (2*m_eff)^(3/2)/2/pi^2/hb^3*sqrt(q*E);
n_mu    = zeros(1, length(mu));

for i = 1 : length(mu)
    n_mu(i) = sum(dos_e./(1+exp(q*(E-mu(i))/(kB*T))))*(dE*q)*1E-6;
end

[~, ind] = min(log(abs(n_mu-n_dop)));
mu       = mu(ind);

end
