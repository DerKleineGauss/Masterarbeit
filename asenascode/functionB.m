function [B] = functionB(r,q,a, U, Lq, g, w, W0, n, delta, L_D)

B= (Potential(r+0.5*q,a,U, g, w, L_D)...
    -Potential(r-0.5*q,a,U, g, w, L_D));%...
    %+1i*functionW(q, Lq/2-delta, W0, n, delta);


end

