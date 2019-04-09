function [Nv, VX, VY, K, EToV, BCType] = rectangularGrid(L_q, L_r)
Globals2D;
% q -> x

h_x = L_q / Kx;
h_y = L_r / Ky;

[VX,VY] = meshgrid(-L_q/2:h_x:L_q/2, -L_r/2:h_y:L_r/2);
VX = VX(:) ; VY = VY(:);

K  = Kx*Ky;
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

Nv = length(VX);

% create BCType
BCType = 0*EToV; % default is 0
tol = 1e-12;
mapSides =  find(abs(abs(VX)-max(VX(:))) < tol);
mapTopBot = find(abs(abs(VY)-max(VY(:))) < tol & abs(abs(VX)-max(VX(:))) > tol);
CorrectBCTable(mapSides, Dirichlet);
CorrectBCTable(mapTopBot, Dirichlet);

VX = VX'; VY = VY';
%%%% plot grid
% trimesh(EToV,VX,VY,zeros(Nv), 'EdgeColor', 'b', 'Linewidth', 2 )
% title ( sprintf ( 'Grid' ) );
% view(2), axis equal, axis off, drawnow
%%%% plot boundary nodes in red
plot(VX,VY,'k.',VX(mapSides), VY(mapSides), 'ro', VX(mapTopBot), VY(mapTopBot), 'bo' );

return