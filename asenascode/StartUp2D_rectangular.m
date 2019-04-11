% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
if (Nx ~= Ny)
    error('Currently, no support for different polynomial degrees in x and y direction');
end
Npx = Nx+1; Npy = Ny+1;
Nfp = Npx;  % Nfpx = Nfpy required (watch error message above)
Np = Npx*Npy; Nfaces=4; NODETOL = 1e-12;

% Compute nodal set
[r,s] = Nodes2D_rectangular(Nx, Ny);

% figure(2)
% plot(r,s,'*')
% grid on

% Build reference element matrices
V = Vandermonde2D_rectangular(Nx, Ny,r,s); invV = inv(V);
MassMatrix = invV'*invV;
[Dr,Ds] = Dmatrices2D_rectangular(Npx, Npy, r,s,V);

% build coordinates of all the nodes
va = EToV(:,1)';
vc = EToV(:,2)';
vb = EToV(:,4)';  % change here!
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

% find all the nodes that lie on each edge
fmask1   = find( abs(s+1) < NODETOL)';
fmask2   = find( abs(r-1) < NODETOL)';  % change here
fmask3   = find( abs(s-1) < NODETOL)';  % change here
fmask4   = find( abs(r+1) < NODETOL)';  % change here
Fmask  = [fmask1;fmask2;fmask3;fmask4]';
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);

% Create surface integral terms
LIFT = Lift2D_rectangular();

% calculate geometric factors
[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);    % J erfolgreich getestet als Flächenverhältnis

% calculate geometric factors
[nx, ny, sJ] = Normals2D_rectangular();
Fscale = sJ./(J(Fmask,:));

% Build connectivity matrix
[EToE, EToF] = tiConnect2D_rectangular(EToV);

% Build connectivity maps
BuildMaps2D;

% Compute weak operators (could be done in preprocessing to save time)
[Vr, Vs] = GradVandermonde2D(N, r, s);
Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V');
