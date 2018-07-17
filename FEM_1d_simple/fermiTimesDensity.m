function [f] = fermiTimesDensity(E, mu)
    global beta E_c;
    f = sqrt(E-E_c) ./ (1+exp(beta * (E-mu))) ;
end