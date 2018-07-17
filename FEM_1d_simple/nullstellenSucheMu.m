function [nullstellen] = nullstellenSucheMu(mu)
global N_D;
nullstellen = N_D - n_electron(mu);
end