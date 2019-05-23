function [V2Dr,V2Ds] = GradVandermonde2D_rectangular(params)

% function [V2Dr,V2Ds] = GradVandermonde2D_rectangular(N,r,s)
% Purpose : Initialize the gradient of the modal basis (i,j) at (r,s) at order Nx, Ny

Npx = params.Npx;
Npy = params.Npy;
r = params.r;
s = params.s;
Np = params.Np;
V2D = zeros(Np,Np);

V2Dr = zeros(length(r),Np);
V2Ds = zeros(length(r),Np);

% Initialize matrices
sk = 1;
for i=0:Npx-1
  for j=0:Npy-1
    [V2Dr(:,sk),V2Ds(:,sk)] = GradSimplex2DP_rectangular(r,s,i,j);
    sk = sk+1;
  end
end
return;
