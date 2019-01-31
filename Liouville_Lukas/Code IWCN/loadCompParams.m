% COMPUTATIONAL DOMAIN PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 66E-9;
Nx = 2*66+1;
Nc = Nx-1;
x  = linspace(-Lx/2, +Lx/2, Nx);
dx = x(2)-x(1);
save('x.mat', 'x');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINITE VOLUME 

Ly = 130E-9;
Ny = 251;
y  = linspace(-Ly/2, +Ly/2, Ny);

y_cell_centers = (y(1:Ny-1)+y(2:Ny))/2;

save('y.mat', 'y');

dy = y(2)-y(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPLEX ABSORBING POTENTIAL IN THE BOUNDARY REGIONS         

W_0              = +1.0;
w_L              = ceil(0.1*Ny);

W                = zeros(1, (Ny-1)/2);
W((Ny-1)/2-w_L:(Ny-1)/2) = W_0*(((0:w_L)).^2./(w_L).^2);
W                = fliplr(W);
W                = [W,0,fliplr(W)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%