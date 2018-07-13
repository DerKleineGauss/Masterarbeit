k_max=.5e10;
k=linspace(0,k_max,5000);
e = 1.60217662e-19;
mu=0.43*e;
m=9.1e-31;
kb=1.38064852e-23;
T=300;
hbar=6.626070040e-34/(2*pi);
y = fermi_dirac(k,mu,T,m,kb,hbar);


f_hat = fft(y);
n = size(f_hat,2);
lin = 1:n;
y_back = zeros(n,1);
for j=1:n
    for k_index=1:n
        y_back(j) = y_back(j) + 1/n * f_hat(k_index)*exp(2*pi*1i*(j-1)*(k_index-1)/n);
    end
end
plot(k,y);
hold on
plot(k,y_back);
hold off


q = 6e-8;
%for q=logspace(20,23)
    z = cos(q*k).*fermi_dirac(k,mu,T,m,kb,hbar);
    %z = fermi_dirac(x,mu,T,m,kb,hbar);
    plot(k,z);
%end
f_hut = 2/(2*pi)*integral(@(value) fermi_dirac_ft(value, q, mu,T, m, kb, hbar),0,1e10)
for tesla=logspace(1,15)
    if fermi_dirac(tesla,mu, T, m, kb, hbar) < 0.99
        k_min=tesla
        break
    end
end
f_hut = 2/(2*pi)* ( integral(@(value) fermi_dirac_ft(value, q, mu,T, m, kb, hbar),k_min,1e10)  )
k2=linspace(k_min,k_max,1000);
plot(k2,fermi_dirac_ft(k2, q, mu,T, m, kb, hbar));

%semilogx(x,y);