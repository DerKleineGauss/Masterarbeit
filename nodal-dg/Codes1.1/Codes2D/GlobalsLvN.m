% Purpose: declare specialized global variables

global m kb e hbar Temperature beta_ N_D E_c;
global L_q L_r

L_q = 152e-9/6;       % nm, lukas paper
L_r = 106e-9/2;       % nm, lukas paper
m=9.1e-31 *0.063;              % Elektronenmasse
kb=1.38064852e-23;
Temperature=300;
beta_ = 1/(kb*Temperature);
hbar=6.626070040e-34/(2*pi);
e = 1.60217662e-19;     % e-charge
E_c = 1.424*e;          % Leitungsbandkante von GaAs (Valenzbandkante auf 0 gesetzt)
N_D = 1e24;             % Donatorkonzentration
