function [x_fixPunkt] = newtonRaphson(f ,x, params)
    h = 1e-15*x;
    f_strich = (f(x+h, params)-f(x-h, params)) / (2*h);
    x_fixPunkt = x - f(x, params)/f_strich;
    while abs((x_fixPunkt-x)/x) > 1e-15
        x = x_fixPunkt;
        h = 1e-15*x;
        f_strich = (f(x+h, params)-f(x-h, params)) / (2*h);
        x_fixPunkt = x - f(x, params)/f_strich;
    end
end