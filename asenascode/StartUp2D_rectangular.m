function [out_params] = StartUp2D_rectangular(params)
% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
if (params.Nx ~= params.Ny)
    error('Currently, no support for different polynomial degrees in x and y direction');
end
Npx = params.Npx;
Npy = params.Npy;
Nfp = params.Nfp;  % Nfpx = Nfpy required (watch error message above)
Np = params.Np;
Nfaces=params.Nfaces;
NODETOL = params.NODETOL;

% Compute nodal set
[params.r,params.s] = Nodes2D_rectangular(params);

% figure(2)
% plot(r,s,'*')
% grid on

% Build reference element matrices
params.V = Vandermonde2D_rectangular(params); params.invV = inv(params.V);
params.MassMatrix = params.invV'*params.invV;
[params.Dr,params.Ds] = Dmatrices2D_rectangular(params);

% build coordinates of all the nodes
r = params.r;
s = params.s;
va = params.EToV(:,1)';
vb = params.EToV(:,2)';
vc = params.EToV(:,4)';  % change here!
params.x = 0.5*(-(r+s)*params.VX(va)+(1+r)*params.VX(vb)+(1+s)*params.VX(vc));
params.y = 0.5*(-(r+s)*params.VY(va)+(1+r)*params.VY(vb)+(1+s)*params.VY(vc));

% find all the nodes that lie on each edge
fmask1   = find( abs(s+1) < NODETOL)';
fmask2   = find( abs(r-1) < NODETOL)';  % changed here
fmask3   = find( abs(s-1) < NODETOL)';  % changed here
fmask4   = find( abs(r+1) < NODETOL)';  % changed here
params.Fmask  = [fmask1;fmask2;fmask3;fmask4]';
params.Fx = params.x(params.Fmask(:), :);
params.Fy = params.y(params.Fmask(:), :);

% Create surface integral terms
params.LIFT = Lift2D_rectangular(params);

% calculate geometric factors
[params.rx,params.sx,params.ry,params.sy,params.J] = GeometricFactors2D(params);    % J erfolgreich getestet als Fl�chenverh�ltnis

% calculate geometric factors
[params.nx, params.ny, params.sJ] = Normals2D_rectangular(params);
params.Fscale = params.sJ./(params.J(params.Fmask,:));

% Build connectivity matrix
[params.EToE, params.EToF] = tiConnect2D_rectangular(params);

% Build connectivity maps
[params.mapM, params.mapP, params.vmapM, params.vmapP, params.vmapB, params.mapB] = BuildMaps2D(params);

% Compute weak operators (could be done in preprocessing to save time)
[Vr, Vs] = GradVandermonde2D_rectangular(params);
params.Drw = (params.V*Vr')/(params.V*params.V'); params.Dsw = (params.V*Vs')/(params.V*params.V');

out_params = params;
end
