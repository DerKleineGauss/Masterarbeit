function [LIFT] = Lift2D_rectangular(params)

% function [LIFT] = Lift2D_rectangular()
% Purpose  : Compute surface to volume lift term for DG formulation

Np = params.Np;
Nfaces = params.Nfaces;
Nfp = params.Nfp;
Nx = params.Nx;
Ny = params.Ny;
Fmask = params.Fmask;
r = params.r;
s = params.s;
V = params.V;

Emat = zeros(Np, Nfaces*Nfp);

% face 1
faceR = r(Fmask(:,1));
V1D = Vandermonde1D(Nx, faceR);
massEdge1 = inv(V1D*V1D');
Emat(Fmask(:,1),1:Nfp) = massEdge1;

% face 2
faceS = s(Fmask(:,2));
V1D = Vandermonde1D(Ny, faceS);
massEdge2 = inv(V1D*V1D');
Emat(Fmask(:,2),Nfp+1:2*Nfp) = massEdge2;

% face 3
faceR = r(Fmask(:,3));
V1D = Vandermonde1D(Nx, faceR);
massEdge3 = inv(V1D*V1D');
Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = massEdge3;

% face 3
faceS = s(Fmask(:,4));
V1D = Vandermonde1D(Ny, faceS);
massEdge4 = inv(V1D*V1D');
Emat(Fmask(:,4),3*Nfp+1:4*Nfp) = massEdge4;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = V*(V'*Emat);
return
