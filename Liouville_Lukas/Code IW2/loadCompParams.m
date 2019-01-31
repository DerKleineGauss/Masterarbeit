% COMPUTATIONAL DOMAIN PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 66E-9;
Nx = 132;
Nc = Nx-1;
x  = linspace(-Lx/2, +Lx/2, Nx);
dx = x(2)-x(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINITE VOLUME 

Ly = 260E-9;
Ny = 134;
y  = linspace(-Ly/2, +Ly/2, Ny);

dy = y(2)-y(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPLEX ABSORBING POTENTIAL IN THE BOUNDARY REGIONS         

W_0              = -1.0;
w_L              = ceil(0.1*Ny);

W                = zeros(1, Ny/2);
W(Ny/2-w_L:Ny/2) = W_0*(((0:w_L)).^2./(w_L).^2);
W                = fliplr(W);
W                = fliplr(W);
Nk               = Ny;

Wr               = zeros(Nk, Nk);
for m = 1 : Nk
    Wr(m,1:Nk) = +q/hb*2/Ny*W*cos(2*pi/Ny*(1:Ny/2).'*((1:Nk)-m));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%