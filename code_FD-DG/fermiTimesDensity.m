function [f] = fermiTimesDensity(E, mu, params)
% E_c = params.E_c;
E_c = 0;
beta = 1/params.constants.kB/params.constants.Temp;
    f = sqrt(E-E_c) ./ (1+exp(beta * (E-mu))) ;
end
