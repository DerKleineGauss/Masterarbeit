function [Nv, VX, VY, K, EToV, BCType] = regularGrid(h)
Globals2D;
GlobalsLvN;
% r -> x
N_x = floor(L_r / h)+1;
h_x = L_r / N_x;
N_y = floor(L_q / h_x)+1;
h_y = L_q / N_y;

[VX,VY] = meshgrid(-L_r/2:h_x:L_r/2, -L_q/2:h_y:L_q/2);
VX = VX(:) ; VY = VY(:);
EToV = delaunay(VX,VY);      
K  = size(EToV,1);  Nv = length(VX);

% triplot(EToV,VX,VY, 'Linewidth', 0.5 )
% title ( sprintf ( 'Grid' ) );

% Reorder elements to ensure counter clockwise orientation
ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
bx = VX(EToV(:,2)); by = VY(EToV(:,2));
cx = VX(EToV(:,3)); cy = VY(EToV(:,3));

D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
i = find(D<0);
EToV(i,:) = EToV(i,[1 3 2]);

% create BCType
BCType = 0*EToV; % default is 0
tol = 1e-12;
mapSides =  find(abs(abs(VX)-max(VX(:))) < tol);
mapTopBot = find(abs(abs(VY)-max(VY(:))) < tol & abs(abs(VX)-max(VX(:))) > tol);
BCType = CorrectBCTable(EToV, BCType, mapSides, Dirichlet);
BCType = CorrectBCTable(EToV, BCType, mapTopBot, Dirichlet);

VX = VX'; VY = VY';
%%%% plot grid
% trimesh(EToV,VX,VY,zeros(Nv), 'EdgeColor', 'b', 'Linewidth', 2 )
% title ( sprintf ( 'Grid' ) );
% view(2), axis equal, axis off, drawnow
%%%% plot boundary nodes in red
plot(VX,VY,'k.',VX(mapSides), VY(mapSides), 'ro', VX(mapTopBot), VY(mapTopBot), 'bo' );

return