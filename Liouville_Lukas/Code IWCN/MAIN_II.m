%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FILE FOR THE NUMERICAL SOLUTION OF THE LIOUVILLE VON NEUMANN EQUA- %
% TION INCLUDING AS WELL DIFFERENT SETS OF BASIS FUNCTIONS                %
% - STATIONARY FLATBAND CASE                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION OF MATLAB
clear all, close all, clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER FILE INCLUDING THE COMPUTATIONAL DOMAIN AND THE CONSTANTS
loadConstants;
loadCompParams;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET THE BASIS
Nk      = ceil(Ny-1);
k       = -2*pi/Ly*(-Nk/2+1/2:+Nk/2-1/2);
dk      = abs(k(2)-k(1));

for i = 1 : Nk
    P(1:Ny-1,i) = exp(1i*k(i)*y_cell_centers)*sqrt(dy/Ly);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFUSION OPERATOR IN THE REAL SAPCE REPRESENTATION
A = spdiags(ones(Ny-1,1)*[-1,+1], [-1,+1], Ny-1, Ny-1)/2;
% EIGENVALUE EXPANSION TECHNIQUE
% [P, lmda] = eig(full(A), 'vector');
% [val, idx] = sort(imag(lmda), 'descend');
% lmda = 1i*val; P = P(:,idx);

% SPAN OF A SUBSPACE
num_of_modes = 60;
P = P(:, Nk/2-num_of_modes/2+1:Nk/2+num_of_modes/2);
k = k(Nk/2-num_of_modes/2+1:Nk/2+num_of_modes/2);
Nk = num_of_modes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AND THE CORRESPONDING BOUNDARY CONDITIONS

fermi_in_k_space = (4*pi*me*kB*T/h^2*log(1+exp(-1/(kB*T)*(h^2*(k).^2/(8*pi^2*me)- muf)))).';
fermi_in_r_space = dk/2/pi*cos(kron(k,y_cell_centers'))*fermi_in_k_space;
% fermi            = P\fermi_in_r_space;
fermi = P'*fermi_in_r_space;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESONANT TUNNELING DEVICE
V_0                    = zeros(1, Nx);
V_0(abs(x) < 7E-9)     = (0.2098+4.7218e-005 +1.2009e-011);%(Eg_AlGaAs-Eg_GaAs)*rel;
V_0(abs(x) < 3E-9)     = 0.0;

V_bias                 = zeros(1, Nx);
V_bias(abs(x) < 13E-9) = linspace(0,V_apl, sum(abs(x) < 13E-9));
V_bias(x >= 13E-9)     = V_apl;

V_stat                      = V_0;
V_tran                      = V_0+V_bias;
% THE CORRESPONDING DOPING PROFIL

Nd_V                 = ones(1, Nx)*Nd;
Nd_V(abs(x) < 13E-9) = 0;

% REFRACTIVE INDEX DISTRIBUTION

eps_V                 = ones(1, Nx)*eps_r_GaAs*eps_0;
eps_V(abs(x) < 7E-9)  = eps_r_AlGaAs*eps_0;
eps_V(abs(x) < 3E-9)  = eps_r_GaAs*eps_0;

eps_V_he  = [eps_V(1), eps_V, eps_V(Nx)];
eps_m     = (eps_V_he(1:Nx)+eps_V_he(2:Nx+1))/2;
eps_p     = (eps_V_he(2:Nx+1)+eps_V_he(3:Nx+2))/2;

% Laplacian = -1/dx/dx*spdiags([eps_p.', -(eps_m+eps_p).', eps_p.'], [-1,0,+1], Nx, Nx);
Laplacian = -1/dx/dx*(diag(eps_p(1:Nx-1).',-1)+diag(-(eps_m+eps_p).',0)+diag(eps_p(1:Nx-1).',+1));
Lap_b     = -1/dx/dx*[eps_V_he(1)*0; zeros(Nx-2,1);eps_V_he(end)*V_apl];
U1 = V_bias'; UU = U1; temp_ch = 1;
while temp_ch > 1E-3
    U1 = U1 + 0.3*(UU-U1);
    V_stat = V_0' - U1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROTATED DEVICE DOMAIN
B_stat = zeros(Ny, Nx);
B_tran = zeros(Ny, Nx);

for i = 1 : (Ny-1)/2
    for j = 1 : Nx
        B_stat(i,j)          = (getSpatTrafPot(x(j)+1/2*y(i), x, V_stat)-getSpatTrafPot(x(j)-1/2*y(i), x, V_stat));
        B_tran(i,j)          = (getSpatTrafPot(x(j)+1/2*y(i), x, V_tran)-getSpatTrafPot(x(j)-1/2*y(i), x, V_tran));
    end
end

B_stat((Ny-1)/2+2:Ny,:) = -flipud(B_stat(1:(Ny-1)/2,:));
B_tran((Ny-1)/2+2:Ny,:) = -flipud(B_tran(1:(Ny-1)/2,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SYSTEM MATRIX ASSEMBLING

L      = sparse((Nx-1)*Nk, Nx*Nk);
R_stat = sparse((Nx-1)*Nk, Nx*Nk);
R_tran = sparse((Nx-1)*Nk, Nx*Nk);

At = +1i*hb/me*P'*A*P;
% At = +1i*hb/me*spdiags(lmda, 0, Nk, Nk);
It = P'*P;
B_stat  = B_stat - 1i*kron(ones(1,size(B_stat,2)),W');
B_tran  = B_tran - 1i*kron(ones(1,size(B_tran,2)),W');

for i = 1 : Nx-1
    Bt1_stat = -1i*q/hb*(P'*spdiags([B_stat(1:Ny-1,i),B_stat(1:Ny-1,i)+B_stat(2:Ny,i),B_stat(2:Ny,i)],-1:1, Ny-1,Ny-1)*dy/4)*P;
    Bt2_stat = -1i*q/hb*(P'*spdiags([B_stat(1:Ny-1,i+1),B_stat(1:Ny-1,i+1)+B_stat(2:Ny,i+1),B_stat(2:Ny,i+1)],-1:1, Ny-1,Ny-1)*dy/4)*P;
    
    Bt1_tran = -1i*q/hb*(P'*spdiags([B_tran(1:Ny-1,i),B_tran(1:Ny-1,i)+B_tran(2:Ny,i),B_tran(2:Ny,i)],-1:1, Ny-1,Ny-1)*dy/4)*P;
    Bt2_tran = -1i*q/hb*(P'*spdiags([B_tran(1:Ny-1,i+1),B_tran(1:Ny-1,i+1)+B_tran(2:Ny,i+1),B_tran(2:Ny,i+1)],-1:1, Ny-1,Ny-1)*dy/4)*P;
    %     Bt1 = -1i*q/hb*(P\spdiags([B(1:Ny-1,i),B(1:Ny-1,i)+B(2:Ny,i),B(2:Ny,i)],-1:1, Nk,Nk)*dy/4)*P;
    %     Bt2 = -1i*q/hb*(P\spdiags([B(1:Ny-1,i+1),B(1:Ny-1,i+1)+B(2:Ny,i+1),B(2:Ny,i+1)],-1:1, Nk,Nk)*dy/4)*P;
    
    
    % LEFT HAND SIDE OF THE SYSTEM MATRIX
    L((i-1)*Nk+1:i*Nk, (i-1)*Nk+1:i*Nk) = It/2*dy;
    L((i-1)*Nk+1:i*Nk, i*Nk+1:(i+1)*Nk) = It/2*dy;
    
    %     Ba = (Bt1+Bt2)/2;
    % RIGHT HAND SIDE OF THE SYSTEM MATRIX
    R_stat((i-1)*Nk+1:i*Nk, (i-1)*Nk+1:i*Nk) = -At/dx+Bt1_stat/2;
    R_stat((i-1)*Nk+1:i*Nk, i*Nk+1:(i+1)*Nk) = +At/dx+Bt2_stat/2;
    
    R_tran((i-1)*Nk+1:i*Nk, (i-1)*Nk+1:i*Nk) = -At/dx+Bt1_tran/2;
    R_tran((i-1)*Nk+1:i*Nk, i*Nk+1:(i+1)*Nk) = +At/dx+Bt2_tran/2;
    %     R((i-1)*Nk+1:i*Nk, (i-1)*Nk+1:i*Nk) = -At/dx+Ba/2;
    %     R((i-1)*Nk+1:i*Nk, i*Nk+1:(i+1)*Nk) = +At/dx+Ba/2;
    
end


bc_r    = R_stat(:, 1:Nk/2);
bc_l    = R_stat(:, (Nx-1)*Nk+Nk/2+1:Nk*Nx);
Rs_red  = R_stat(:, Nk/2+1:(Nx-1)*Nk+Nk/2);clear R_stat
Rt_red  = R_tran(:, Nk/2+1:(Nx-1)*Nk+Nk/2);clear R_tran
L_red   = L(:, Nk/2+1:(Nx-1)*Nk+Nk/2);clear L
bc      = (bc_r*fermi(1:Nk/2)+bc_l*fermi(Nk/2+1:Nk));clear bc_r, clear bc_l

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE STATIONARY SOLUTION FOR THE EXPANSION COEFFICIENT

a    = (Rs_red)\(-bc);clear Rs_Red
af   = [fermi(1:Nk/2); a; fermi(Nk/2+1:Nk)];
af   = full(reshape(af, Nk, Nx));
C    = P*af;

figure(1), mesh(real(C)); title(num2str(i)); drawnow;
figure(2), mesh(real(af)); title(num2str(i)); drawnow;
figure(3), mesh(imag(C)); title(num2str(i)); drawnow;

n = real(C((Ny-1)/2,:)+C((Ny-1)/2+1,:))/2;
figure(2), plot(n), drawnow;

%  UU=-Laplacian\Lap_b+Laplacian\(n'-Nd_V')*q;
 UU = getHartreeFockPotential((eps_V_he(1:end-1)+eps_V_he(2:end))/2, U1, n, Nx, dx, Nd_V);
 temp_ch=max(max((abs(UU-U1))));
 figure(41), plot(V_stat), drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
