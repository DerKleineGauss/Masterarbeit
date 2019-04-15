function [V2D] = Vandermonde2D(Nx, Ny, r, s);

% function [V2D] = Vandermonde2D(N, r, s);
% Purpose : Initialize the 2D Vandermonde Matrix,  V_{ij} = phi_j(r_i, s_i);

V2D = zeros(length(r), Nx*Ny);


% build the Vandermonde matrix
sk = 1;
for i=0:Nx
  for j=0:Ny
    V2D(:,sk) = Simplex2DP(r,s,i,j);
    sk = sk+1;
  end
end
return;
