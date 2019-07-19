classdef PhysicalConstants
  properties (Constant)
    Lr= 106e-9/2;
    Lq= 152e-9/2;
    L_D= 40e-9/2;
    w= 4e-9;
    g= 6e-9;
    hbar= 1.054571800e-34;
    e= 1.6021766208e-19;
    W0= 1;
    n= 1;
    a0= 0.1768 ;
    m0= 9.10938356e-31;
    kB= 1.38064852e-23;
    Temp= 300;
    N_D = 1e24;             % Donatorkonzentration
  end
  properties
    E_c = 0;                % Leitungsbandkantenenergie
    delta= 0.2*PhysicalConstants.Lq;
    m= 0.063*PhysicalConstants.m0;
    mu_l= 0.0467*PhysicalConstants.e;
    mu_r= 0.0467*PhysicalConstants.e;
  end
end
