x           = 0.33;
m0          = 9.109E-31;
q           = 1.602E-19;
rel         = 0.57;
h           = 6.626E-34;
hb          = h/2/pi;
kB          = 1.38E-23;
T           = 300;
me          = 0.063*m0;
Nd          = 1E+24;
Nd2         = 0;
muf         = q*getFermiLevel(Nd/1E+6, me, T);
Eg_AlGaAs   = 1.424+1.721*x-1.437*x^2+1.31*x^3;
Eg_GaAs     = 1.424;

eps_0        =   8.8541878176e-12; 
eps_r_GaAs   =   12.9; 
eps_r_AlGaAs =   (12.9-2.84*x); 
