function [nx, ny, sJ] = Normals2D_rectangular()

% function [nx, ny, sJ] = Normals2D_rectangular()
% Purpose : Compute outward pointing normals at elements faces and surface Jacobians

Globals2D;
xr = Dr*x; yr = Dr*y; xs = Ds*x; ys = Ds*y; J = xr.*ys-xs.*yr;

% interpolate geometric factors to face nodes
fxr = xr(Fmask, :); fxs = xs(Fmask, :); fyr = yr(Fmask, :); fys = ys(Fmask, :);

% build normals
nx = zeros(3*Nfp, K); ny = zeros(3*Nfp, K);
fid1 = (1:Nfp)'; fid2 = (Nfp+1:2*Nfp)'; fid3 = (2*Nfp+1:3*Nfp)'; fid4 = (3*Nfp+1:4*Nfp)';

% face 1
nx(fid1, :) =  fyr(fid1, :); ny(fid1, :) = -fxr(fid1, :);

% face 2
nx(fid2, :) =  fys(fid2, :); ny(fid2, :) = -fxs(fid2, :);

% face 3
nx(fid3, :) = -fyr(fid3, :); ny(fid3, :) =  fxr(fid3, :);

% face 4
nx(fid4, :) = -fys(fid4, :); ny(fid4, :) =  fxs(fid4, :);

nx(abs(nx/max(nx(:))) < 1e-14)=0;
ny(abs(ny/max(ny(:))) < 1e-14)=0;

% normalise
sJ = sqrt(nx.*nx+ny.*ny); nx = nx./sJ; ny = ny./sJ;
return;
