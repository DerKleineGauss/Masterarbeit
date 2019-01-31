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
num_of_modes = Nk;
P = P(:, Nk/2-num_of_modes/2+1:Nk/2+num_of_modes/2);
k = k(Nk/2-num_of_modes/2+1:Nk/2+num_of_modes/2);
Nk = num_of_modes;

save('k.mat', 'k');

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
save('V_stat.mat','V_stat');

Nd_V                 = ones(1, Nx)*Nd;
Nd_V(abs(x) < 13E-9) = 0;

save('N.mat','Nd_V');

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

a    = (Rt_red)\(-bc);clear Rs_Red
af   = [fermi(1:Nk/2); a; fermi(Nk/2+1:Nk)];
af   = full(reshape(af, Nk, Nx));
C    = P*af;

figure(1), mesh(real(C)); title(num2str(i)); drawnow;
figure(2), mesh(real(af)); title(num2str(i)); drawnow;
figure(3), mesh(imag(C)); title(num2str(i)); drawnow;
n    = real(C((Ny-1)/2,:)+C((Ny-1)/2+1,:))/2;
save(['C', num2str(num_of_modes), '.mat'], 'C');
save(['n', num2str(num_of_modes),'.mat'], 'n');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % GET THE TRANSIENT SOLUTION FOR THE EXPANSION COEFFICIENT REGARDING THE
% % % NONEQUILIBRIUM SOLUTION
% % 
% % % Nt = 801;
% % % t  = linspace(0, 0.25E-13, Nt);
% % % dt = t(2)-t(1);
% % % 
% % % dt = 3.125E-17;
% % 
% % % CFL-CONDITION;
% % 
% % dt = 0.05*(dx/max(hb*k)*me);
% % % dt = 1E-17;
% % Nt = 20001;
% % 
% % U  = L_red\Rt_red;
% % b0 = L_red\bc;
% % 
% % for i = 1 : Nt
% %     
% %     a    = getLSRK_4(U,a,dt,b0);
% %     af   = [fermi(1:Nk/2); a; fermi(Nk/2+1:Nk)];
% %     af   = full(reshape(af, Nk, Nx));
% %     C    = P*af;
% %     figure(1), mesh(real(C)); title(num2str(i)); drawnow;
% %     figure(2), imagesc(real(af)); title(num2str(i)); drawnow;
% % %     figure(3), mesh(imag(C)); title(num2str(i)); drawnow;
% %     
% % end
% % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %
% % % a    = [fermi(1:Nk/2); a; fermi(Nk/2+1:Nk)];
% % % af   = full(reshape(a, Nk, Nx));
% % % % n    = abs(dk)/2/pi*sum(af,1);
% % % % j_o  = -q*abs(dk)/2/pi*hb/me*(k*af(:,1:Nx));
% % %
% % % j_o = -q*hb/me*(k*af(:,1:Nx))*sqrt(dy/Ly);
% % %
% % % C    = P*af;
% % %
% % % figure, mesh(real(C));
% % % figure, mesh(real(af));
% % % figure, mesh(imag(C));
% % %
% % % J_op = 1i*q*hb/me*kron(speye(Nx,Nx), spdiags([-1/2/dy*ones(Ny-1,1),+1/2/dy*ones(Ny-1,1)],[-1, +1], Ny-1,Ny-1));
% % % j = J_op*C(:);
% % % j = reshape(j, Ny-1, Nx);
% % %
% % % figure, mesh(real(j));
% % %
% % % j = real(j((Ny-1)/2,:)+j((Ny-1)/2+1,:))/2;
% % % j = -q*hb/me*imag(+C((Ny-1)/2+1,:)-C((Ny-1)/2,:))/dy;
% % % figure, plot(real(j_o), 'color', 'r'), hold on, plot(j, 'color', 'b');
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
