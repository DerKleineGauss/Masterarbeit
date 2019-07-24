classdef GlobalParams
  properties
    % Purpose: declare global variables
    constants = PhysicalConstants;
    makeMovie
    movieName
    saveResults
    plot_logarithmic
    mode
    timestepping
    withCAP
    U   % voltage
    Np
    Npy % number of interfaces coordinates in y space originating from FD method
    Ny  % number of cells in y space originating from FD method
    Nfp
    N
    nx  % normals in x
    Nfaces
    K
    r
    rx
    Dr
    LIFT
    MassMatrix
    Fx
    Fscale
    J
    Fmask
    loc2glb
    mapB
    mapI
    mapO
    mapW
    mapF
    mapC
    mapS
    mapM
    mapP
    mapD
    mapN
    vmapM
    vmapP
    vmapB
    vmapI
    vmapO
    vmapW
    vmapC
    vmapS
    vmapD
    vmapN
    vmapF
    rk4a = [            0.0 ...
            -567301805773.0/1357537059087.0 ...
            -2404267990393.0/2016746695238.0 ...
            -3550918686646.0/2091501179385.0  ...
            -1275806237668.0/842570457699.0];
    rk4b = [ 1432997174477.0/9575080441755.0 ...
             5161836677717.0/13612068292357.0 ...
             1720146321549.0/2090206949498.0  ...
             3134564353537.0/4481467310338.0  ...
             2277821191437.0/14882151754819.0];
    rk4c = [             0.0  ...
             1432997174477.0/9575080441755.0 ...
             2526269341429.0/6820363962896.0 ...
             2006345519317.0/3224310063776.0 ...
             2802321613138.0/2924317926251.0];
    EToE
    EToF
    EToV
    V
    invV
    x
    x_interface
    y
    y_interface
    hy
    hx
    NODETOL = 1e-12;
    VX
    testing
    FinalTime
    rampTime
    % Some curved mesh and cubature quadrature specific data
    cub
    gauss
    straight
    curved
    % boundary maps
    In = 1;
    Out = 2;
    Wall = 3;
    Far = 4;
    Cyl = 5;
    Dirichlet = 6;
    Neuman = 7;
    Slip = 8;
    % scaled values
    Lr_scaled
    Lq_scaled
    L_D_scaled
    w_scaled
    g_scaled
    delta_scaled
    epsilon
    gamma
  end
  methods
    function [Lr, Lq, L_D, w, g, delta] = scale(obj, gamma)
      Lr= floor(obj.constants.Lr*gamma);
      Lq= floor(obj.constants.Lq*gamma);
      L_D= obj.constants.L_D*gamma;
      w= obj.constants.w*gamma;
      g= obj.constants.g*gamma;
      delta= obj.constants.delta*gamma;
    end
    function [a, b, c_l, c_r] = get_abc(obj, gamma, epsilon)
      a= gamma*obj.constants.m*obj.constants.kB*obj.constants.Temp/2/(pi*obj.constants.hbar)^2;
      b= epsilon*obj.constants.hbar/obj.constants.kB/obj.constants.Temp / 2;
      c_l= obj.constants.mu_l/obj.constants.kB/obj.constants.Temp;
      c_r= obj.constants.mu_r/obj.constants.kB/obj.constants.Temp;
    end
    function [t_fs] = characteristicTimeToFs(obj, t_ch)
      t_fs = t_ch / obj.epsilon * 1e15;
    end
  end
end
