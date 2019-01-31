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
Nk      = Ny;
k       = -2*pi/Ly*(-Nk/2+1/2:+Nk/2-1/2);
dk      = abs(k(2)-k(1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFUSION OPERATOR IN THE REAL SAPCE REPRESENTATION
A = spdiags(k.', 0, Nk, Nk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AND THE CORRESPONDING BOUNDARY CONDITIONS

fermi_in_k_space = (4*pi*me*kB*T/h^2*log(1+exp(-1/(kB*T)*(h^2*(k).^2/(8*pi^2*me)- muf)))).';

fermi            = fermi_in_k_space;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESONANT TUNNELING DEVICE
V_0                    = zeros(1, Nx);
V_0(abs(x) < 7E-9)     = (0.2098+4.7218e-005 +1.2009e-011);%(Eg_AlGaAs-Eg_GaAs)*rel;
V_0(abs(x) < 3E-9)     = 0.0;

V_bias                 = zeros(1, Nx);
V_bias(abs(x) < 13E-9) = linspace(0,V_apl, sum(abs(x) < 13E-9));
V_bias(x >= 13E-9)     = V_apl;

V_stat                 = V_0;
V_tran                 = V_0+V_bias;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROTATED DEVICE DOMAIN
B_stat = zeros(Ny/2, Nx);
B_tran = zeros(Ny/2, Nx);

for i = 1 : Ny/2
    for j = 1 : Nx
        B_stat(i,j)          = (getSpatTrafPot(x(j)+1/2*y(i), x, V_stat)-getSpatTrafPot(x(j)-1/2*y(i), x, V_stat));
        B_tran(i,j)          = (getSpatTrafPot(x(j)+1/2*y(i), x, V_tran)-getSpatTrafPot(x(j)-1/2*y(i), x, V_tran));
    end
end

B_stat = flipud(B_stat);
B_tran = flipud(B_tran);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SYSTEM MATRIX ASSEMBLING

L      = sparse((Nx-1)*Nk, Nx*Nk);
R_stat = sparse((Nx-1)*Nk, Nx*Nk);
R_tran = sparse((Nx-1)*Nk, Nx*Nk);

At = -hb/me*A/dx;
It = speye(Nk, Nk);

for i = 1 : Nx-1
    Br1_stat = sparse(Nk, Nk);
    Br2_stat = sparse(Nk, Nk);
    Br1_tran = sparse(Nk, Nk);
    Br2_tran = sparse(Nk, Nk);
    
    for m = 1 : Nk
        Br1_stat(m,1:Nk) = +q/hb*2/Ny*(B_stat(:,i)).'*sin(2*pi/Ny*(1:Ny/2).'*((1:Nk)-m));
        Br2_stat(m,1:Nk) = +q/hb*2/Ny*(B_stat(:,i+1)).'*sin(2*pi/Ny*(1:Ny/2).'*((1:Nk)-m));
        Br1_tran(m,1:Nk) = +q/hb*2/Ny*(B_tran(:,i)).'*sin(2*pi/Ny*(1:Ny/2).'*((1:Nk)-m));
        Br2_tran(m,1:Nk) = +q/hb*2/Ny*(B_tran(:,i+1)).'*sin(2*pi/Ny*(1:Ny/2).'*((1:Nk)-m));
    end
    
    
    % LEFT HAND SIDE OF THE SYSTEM MATRIX
    L((i-1)*Nk+1:i*Nk, (i-1)*Nk+1:i*Nk) = It/2;
    L((i-1)*Nk+1:i*Nk, i*Nk+1:(i+1)*Nk) = It/2;
    
    %     Ba = (Bt1+Bt2)/2;
    % RIGHT HAND SIDE OF THE SYSTEM MATRIX
    R_stat((i-1)*Nk+1:i*Nk, (i-1)*Nk+1:i*Nk) = -At+(Br1_stat+Wr)/2;
    R_stat((i-1)*Nk+1:i*Nk, i*Nk+1:(i+1)*Nk) = +At+(Br2_stat+Wr)/2;
    
    R_tran((i-1)*Nk+1:i*Nk, (i-1)*Nk+1:i*Nk) = -At+(Br1_tran+Wr)/2;
    R_tran((i-1)*Nk+1:i*Nk, i*Nk+1:(i+1)*Nk) = +At+(Br2_tran+Wr)/2;
    
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
n    = abs(dk)/2/pi*sum(af,1);
j_o  = -q*abs(dk)/2/pi*hb/me*(k*af(:,1:Nx));

%save(['D:\Lukas Schulz\Quantum_Transport_Equations\CONFERENCES\IW2-2019\Code\Code\Transient_Results\af_', num2str(0), 'femtos.mat'], 'af');


a_conv    = (Rt_red)\(-bc);
af_conv   = [fermi(1:Nk/2); a_conv; fermi(Nk/2+1:Nk)];
af_conv   = full(reshape(af_conv, Nk, Nx));
n_conv    = abs(dk)/2/pi*sum(af_conv,1);
j_o_conv  = -q*abs(dk)/2/pi*hb/me*(k*af_conv(:,1:Nx));

%save('D:\Lukas Schulz\Quantum_Transport_Equations\CONFERENCES\IW2-2019\Code\Code\Transient_Results\af_inf_femtos.mat', 'af_conv');

figure(2), mesh(real(af)); title(num2str(i)); drawnow;
figure(3), plot(j_o), drawnow;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GET THE TRANSIENT SOLUTION FOR THE EXPANSION COEFFICIENT REGARDING THE
% % NONEQUILIBRIUM SOLUTION
%
% CFL-CONDITION;

% dt = 0.5*(dx/max(hb*k)*me);
dt = 1E-17;
Nt = 800001;
t  = (0:Nt-1)*dt;
U  = L_red\Rt_red;
b0 = L_red\bc;

for i = 1 : Nt
    
    a    = getLSRK_4(U,a,dt,b0);
    af   = [fermi(1:Nk/2); a; fermi(Nk/2+1:Nk)];
    af   = full(reshape(af, Nk, Nx));
    n    = abs(dk)/2/pi*sum(af,1);
    j_o  = -q*abs(dk)/2/pi*hb/me*(k*af(:,1:Nx));
    
    if mod(i,100)==0
        figure(2), imagesc(real(af)); title(num2str((i)*dt/1E-15)); drawnow;
        figure(3), plot(j_o, 'color','r'), hold on, plot(j_o_conv, 'color','b'), drawnow;
        figure(4), plot(n, 'color', 'r'), hold on, plot(n_conv, 'color', 'b'), drawnow;
        %save(['D:\Lukas Schulz\Quantum_Transport_Equations\CONFERENCES\IW2-2019\Code\Code\Transient_Results\af_', num2str(i*dt/1E-15), 'femtos.mat'], 'af');
    end
end
%
