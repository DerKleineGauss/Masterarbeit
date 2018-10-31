function [x_fixPunkt] = newtonRaphson(f ,x)
    h = 1e-15*x;
    f_strich = (f(x+h)-f(x-h)) / (2*h);
    x_fixPunkt = x - f(x)/f_strich;
    while abs((x_fixPunkt-x)/x) > 1e-15
        x = x_fixPunkt;
        h = 1e-15*x;
        f_strich = (f(x+h)-f(x-h)) / (2*h);
        x_fixPunkt = x - f(x)/f_strich;
    end
end