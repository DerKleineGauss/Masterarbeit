function [V] = Potential(x,a, U, g, w, L_D)

V= a.*(heaviside(x+w/2+g)-heaviside(x-w/2+g)+heaviside(x+w/2-g)-heaviside(x-w/2-g))...
    -U.*heaviside(x-L_D/2).*(-x./L_D+0.5)-U.*heaviside(x+L_D/2).*(x./L_D+0.5);

end

