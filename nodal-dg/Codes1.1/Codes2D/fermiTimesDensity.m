function [f] = fermiTimesDensity(E, mu)
GlobalsLvN;
    f = sqrt(E-E_c) ./ (1+exp(beta_ * (E-mu))) ;
end
