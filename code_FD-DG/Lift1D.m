function [LIFT] = Lift1D(params)

% function [LIFT] = Lift1D
% Purpose  : Compute surface integral term in DG formulation

Emat = zeros(params.Np , params.Nfaces*params.Nfp);

% Define Emat
Emat(1,1) = 1.0; Emat(params.Np,2) = 1.0;

% inv(mass matrix)*\s_n (L_i,L_j)_{edge_n}
LIFT = params.V*(params.V'*Emat);
return
