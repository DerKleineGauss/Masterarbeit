function [nullstellen] = nullstellenSucheMu(mu, params)
N_D = params.constants.N_D;
nullstellen = N_D - n_electron(mu, params);
end
