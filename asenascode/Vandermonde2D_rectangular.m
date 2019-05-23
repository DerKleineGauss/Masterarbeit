function [V2D] = Vandermonde2D_rectangular(params);

% function [V2D] = Vandermonde2D(N, r, s);
% Purpose : Initialize the 2D Vandermonde Matrix,  V_{ij} = phi_j(r_i, s_i);
Nx = params.Nx;
Ny = params.Ny;
r = params.r;
s = params.s;
Np = params.Np;
V2D = zeros(Np,Np);

% build the Vandermonde matrix
sk = 1;
for i=0:Nx
  for j=0:Ny
    V2D(:,sk) = Simplex2DP_rectangular(r,s,i,j);
    sk = sk+1;
  end
end
return;
