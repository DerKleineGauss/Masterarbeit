% CONSTANTS AND MATERIAL PARAMETERS FOR ALGAAS ALLOYS
x_part      = 0.33;
m0          = 9.109E-31;
q           = 1.602E-19;
rel         = 0.57;
h           = 6.626E-34;
hb          = h/2/pi;
kB          = 1.38E-23;
T           = 300;
me          = 0.063*m0;
Nd          = 1E+24;
muf         = q*get_FermiLevel(Nd/1E+6, me, T);

Eg_AlGaAs   = 1.424+1.721*x_part-1.437*x_part^2+1.31*x_part^3;
Eg_GaAs     = 1.424;

eps_0        =   8.854E-12; 
eps_r_GaAs   =   12.9; 
eps_r_AlGaAs =   (12.9-2.84*x_part); 

% EXTERNAL APPLIED FIELDS
V_apl        = +0.05;

