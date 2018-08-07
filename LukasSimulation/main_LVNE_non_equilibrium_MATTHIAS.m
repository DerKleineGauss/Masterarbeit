% Numerical Solution of the Liouville von Neumann Equation                %
% - selfconsistent non equilibrium approach                               %
% - different ways of the current evaluation                              %
% - finite volume scheme                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

% get the constants needed for the simulation

getConstants;

% initialize the computational domain

L_r         = 66E-9;
N_inter_r   = 3*66+1;
N_cell_r    = N_inter_r-1;
r_inter     = linspace(-L_r/2, +L_r/2, N_inter_r);
r_cell      = (r_inter(1:N_cell_r)+r_inter(2:N_cell_r+1))/2;
dr          = diff(r_inter);

L_q         = 2*66E-9;
N_inter_q   = 4*66+1;
N_cell_q    = N_inter_q-1;
q_inter     = linspace(-L_q/2, +L_q/2, N_inter_q);
q_cell      = (q_inter(1:N_cell_q)+q_inter(2:N_cell_q+1))/2;
dq          = diff(q_inter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for the test device / to be out sourced....

% potential distribution

V                         = zeros(1, N_cell_r);
V(abs(r_cell) < 8E-9)     = (Eg_AlGaAs-Eg_GaAs)*rel;%(0.2098+4.7218e-005 +1.2009e-011);%
V(abs(r_cell) < 3E-9)     = 0.0;
V(N_cell_r/2+1:N_cell_r)  = fliplr(V(1:N_cell_r/2));

% doping profile

Nd_V                        = ones(1, N_cell_r)*Nd;
Nd_V(abs(r_cell) < 14E-9)   = 0;
Nd_V(N_cell_r/2+1:N_cell_r) = fliplr(Nd_V(1:N_cell_r/2));

% refractive index distribution for Poisson equation

eps_V                        = ones(1, N_inter_r)*eps_r_GaAs*eps_0;
eps_V(abs(r_inter) < 8E-9)    = eps_r_AlGaAs*eps_0;
eps_V(abs(r_inter) < 3E-9)    = eps_r_GaAs*eps_0;
eps_V(N_cell_r/2+1:N_cell_r) = fliplr(eps_V(1:N_cell_r/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the basis of the diffusion operator

D       = zeros(1, N_cell_q);
Q       = zeros(N_cell_q, N_cell_q);

for k = 1 : N_cell_q
    D(k) = 1i*cos(pi*k/(N_cell_q+1));
    for j = 1 : N_cell_q
        Q(j,k) = (1i)^(j)*sin(k*pi*j/((N_cell_q+1)));
    end
end

% complex absorbing potential

W                               = zeros(N_inter_q, 1);
W_0                             = 5;
delta                           = 0.2*L_q;
n_delta                         = ceil(delta/dq(N_cell_q));
W(N_inter_q-n_delta+1:N_inter_q)= W_0*(1:n_delta).^4/n_delta^4;
W(1:n_delta)                    = flipud(W(N_inter_q-n_delta+1:N_inter_q));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of the boundary conditions / to be outsourced.....

% phase space representation
k                           = 2*pi/L_q*(-N_cell_q/2+1/2:N_cell_q/2-1/2);
dk                          = k(2)-k(1);
fermi_in_phase_space        = 4*pi*me*kB*T/h^2*log(1+exp(-1/(kB*T)*(h^2*(k).^2/(8*pi^2*me)- muf))).';
fermi_in_real_space         = dk/2/pi*cos(q_inter.'*k)*fermi_in_phase_space;
fermi_in_real_space         = (fermi_in_real_space(1:N_cell_q)+fermi_in_real_space(2:N_cell_q+1))/2;

fermi_in_eig_space          = Q\fermi_in_real_space;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non Equilibrium Simulation

V_apl                       = +0.27;

V_bias                      = zeros(1, N_cell_r);
V_bias(abs(r_cell) < 14E-9) = linspace(0,V_apl, sum(abs(r_cell) < 14E-9));
V_bias(r_cell >= 14E-9)     = V_apl;
V_H                         = V_bias;

B = zeros(N_inter_q, N_cell_r);
F = sparse(N_cell_q*N_cell_r, N_cell_q*N_inter_r);

% potential in the rotated domain

for i = 1 : (N_inter_q-1)/2
    for j = 1 : N_cell_r
        B(i,j)          = (getPotential_rq(r_cell(j)+1/2*q_inter(i), r_cell, V-V_H)-getPotential_rq(r_cell(j)-1/2*q_inter(i), r_cell, V-V_H));
    end
end
B((N_inter_q-1)/2+1,:)              = 0;
B(((N_inter_q-1)/2+2):N_inter_q,:)  = -flipud(B(1:(N_inter_q-1)/2,:));

B = B-1j*W*ones(1,size(B,2));

for i = 1 : N_cell_r
    G = sparse(N_cell_q, N_cell_q);
    G = (1/4*diag(B(2:N_cell_q, i), -1) + 1/4*diag(B(1:N_cell_q,i)+B(2:N_cell_q+1,i), 0) + 1/4*diag(B(2:N_cell_q, i), +1));
    
    F(N_cell_q*(i-1)+1:N_cell_q*(i), N_cell_q*(i-1)+1:N_cell_q*(i)) = -expm(+me*q/hb^2*diag(1./D)*(Q\G*Q)*dr(1)/2*diag(dq));%-getTaylorApproximation(+me*q/hb^2*diag(1./D)*(Q\G*Q)*dr(1)/2*diag(dq),2);%
    F(N_cell_q*(i-1)+1:N_cell_q*(i), N_cell_q*(i)+1:N_cell_q*(i+1)) = +expm(-me*q/hb^2*diag(1./D)*(Q\G*Q)*dr(1)/2*diag(dq));%+getTaylorApproximation(-me*q/hb^2*diag(1./D)*(Q\G*Q)*dr(1)/2*diag(dq),2);%
    
end

bc_l    = F(:, 1:N_cell_q/2);
bc_r    = F(:, N_cell_r*N_cell_q+N_cell_q/2+1:N_cell_q*N_inter_r);
F       = F(:, N_cell_q/2+1:N_cell_r*N_cell_q+N_cell_q/2);
bc      = (bc_l*fermi_in_eig_space(1:N_cell_q/2)+bc_r*fermi_in_eig_space(N_cell_q/2+1:N_cell_q));
P       = F\(-bc);
P       = reshape([fermi_in_eig_space(1:N_cell_q/2); P; fermi_in_eig_space(N_cell_q/2+1:N_cell_q)], N_cell_q, N_cell_r+1);
PP      = Q*(P(:,1:N_cell_r)+P(:,2:N_cell_r+1))/2;

% print the results

figure(1), mesh(real(PP)), axis tight, drawnow;
figure(2), plot(real(PP(N_cell_q/2,:)+PP(N_cell_q/2+1))/2), hold on, drawnow;
