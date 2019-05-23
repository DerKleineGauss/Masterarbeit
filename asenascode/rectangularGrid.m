function [out_params] = rectangularGrid(params)
% r -> x

Lr = params.Lr_scaled;
Lq = params.Lq_scaled;
K = params.K;
Kx = params.Kx;
Ky = params.Ky;
h_x = Lr/ Kx;
h_y = Lq / Ky;

[VX,VY] = meshgrid(-Lr/2:h_x:Lr/2, -Lq/2:h_y:Lq/2);
VX = VX(:) ; VY = VY(:);

EToV=zeros(K,4);
for ky=1:Ky
    for kx=1:Kx
        k_actual = (ky-1)*Kx + kx;
        EToV(k_actual,1) = kx*(Ky+1) - (ky-1);
    end
end
EToV(:,2) = EToV(:,1) - 1;
EToV(:,3) = EToV(:,2) + (Ky+1);
EToV(:,4) = EToV(:,3) + 1;
params.EToV = EToV;

% create BCType
params.BCType = 0*EToV; % default is 0
tol = 1e-12;
mapSides =  find(abs(abs(VX)-max(VX(:))) < tol);
mapTopBot = find(abs(abs(VY)-max(VY(:))) < tol);
params = CorrectBCTable_rectangular(params, mapSides, params.Dirichlet);
params = CorrectBCTable_rectangular(params, mapTopBot, params.Dirichlet);

params.VX = VX'; params.VY = VY';
%%%% plot grid
% trimesh(params.EToV,params.VX,params.VY,zeros(params.Nv), 'EdgeColor', 'b', 'Linewidth', 2 )
% title ( sprintf ( 'Grid' ) );
% view(2), axis equal, axis off, drawnow
%%%% plot boundary nodes in red
plot(params.VX,params.VY,'k.',params.VX(mapSides), params.VY(mapSides), 'ro', params.VX(mapTopBot), params.VY(mapTopBot), 'bo' );
out_params = params;
return
