function [f] = LvN_f(r,q)
% calculate em/h^2 B(r,q) = em/h^2 (V(r+q/2) - V*(r-q/2))

GlobalsLvN;

V_p = 0*r;
V_m = 0*r;

V_p(abs(r+q/2)>=L1_half & abs(r+q/2)<=(L1_half+L2)) = max_V;
V_m(abs(r-q/2)>=L1_half & abs(r-q/2)<=(L1_half+L2)) = max_V;
% factor 2 because 1/2 div q = f
f = 2* e*m/(hbar^2) * (V_p - V_m) ;

end
