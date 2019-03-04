function Pn = getLSRK_4(A, P, dt, bc)

a = [0, -567301805773/1357537059087, -2404267990393/2016746695238, -3550918686646/2091501179385, -1275806237668/842570457699];
b = [1432997174477/9575080441755, 5161836677717/13612068292357, 1720146321549/2090206949498, 3134564353537/4481467310338, 2277821191437/14882151754819];
% c = [0, 1432997174477/9575080441755, 2526269341429/6820363962896, 2006345519317/3224310063776, 2802321613138/2924317926251];

kn = 1;
for i = 1 : 5
    kn = a(i)*kn + dt*(A*P+bc);
    P = P + b(i)*kn;
end
Pn = P;
end