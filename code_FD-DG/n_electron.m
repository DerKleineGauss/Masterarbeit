function [y] = n_electron(mu, params)

% trapez iteration
% error_ratio = 1e-10;
% % first iteration
% l = 1;
% a = 0;
% b = 1/E_c;
% T_l = (b-a)/2 * (0 + fermiTimesDensity(b, beta, mu, E_c));  %  f�r eps -> 0 bzw. E -> infty dominiert die exp-Funktion im Nenner, der Term wird 0
% T_ll = 1e20;    % some high value to start iteration
% % iterate with increasing l
% while abs(T_ll -T_l)/T_ll > error_ratio
%     if l>1
%         T_l = T_ll;
%     end
%     l = l + 1;
%     h_l = (b-a)/(2^l);
%     h_ll = h_l/2;
%     tmp = 0;
%     for k=0:(2^(l-1) - 1)
%        tmp = tmp + fermiTimesDensity(a+h_ll+k*h_l, beta, mu, E_c);
%     end
%     T_ll = 1/2 * (T_l + h_l*tmp)
% end
% y = T_ll;

% stop integration at E_m with condition:  beta(E_m-mu) > 20 (Nenner
% explodiert)

%% Testing
%test = linspace(E_c, 20/beta+mu, 1000);
%plot(test,fermiTimesDensity(test,mu));
%%
% E_c = params.E_c;
E_c = 0;
beta = 1/params.constants.kB/params.constants.Temp;
m= params.constants.m;
hbar = params.constants.hbar;

y = real(integral( @(E)fermiTimesDensity(E,mu, params), E_c, (20/beta + mu))) * (2*m)^(3/2) / ( 2*pi^2 * hbar^3 );
end
