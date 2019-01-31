% Driver script for solving the 2D Poisson equation
clear()
close all
figurecounter = 1;

Globals2D;
GlobalsLvN;
GlobalsLvN_data;

% Polynomial order used for approximation
N = 4;

%% Read in Mesh , Initialize solver
%h_stepsize = 1;
%[nodes, EToV, BCType, cmdout] = Mesh(h_stepsize, L_r*1e9 ,L_q*1e9, 2*L2*1e9, 2*L1_half*1e9); % alles in nm
%Nv= size(nodes,1);
%K= size(EToV,1);
%VX= 1e-9*nodes(:,1)';
%VY= 1e-9*nodes(:,2)';
% Read in Mesh - Regelm‰ﬂig
h_stepsize = 2;
[Nv, VX, VY, K, EToV, BCType] = regularGrid(h_stepsize*1e-9);
% plotting the rechengebiet
figure(figurecounter)
figurecounter = figurecounter +1;
triplot(EToV,VX,VY, 'Linewidth', 0.2 )
title ( sprintf ( 'Grid' ) );
% plot(x,y,'.');
% hold on
% linecolors = jet(K);
% for i = 1:K
%     fill(x(:,i),y(:,i), linecolors(i,:));
% end
% hold off

% Initialize solver
StartUp2D;
x = x/(L_r/2);
y = y/(L_q/2);

% set up boundary conditions
BuildBCMaps2D;

%% set up right hand side
% evaluate potential term
q_coordinates = y*(L_q/2); % 1d real coordinates
r_coordinates = x*(L_r/2); % 1d real coordinates
pot_vector = LvN_f(r_coordinates,q_coordinates);
% build diagonal matrix representing f(r,q) in the f(r,q)*u part
pot = sparse(1:K*Np, 1:K*Np, pot_vector(:), K*Np, K*Np);
% % plot test
figure(figurecounter)
figurecounter = figurecounter + 1;
[xq,yq] = meshgrid(-1:0.01:1, -1:0.01:1);
% vq = griddata(x,y,pot_vector,xq,yq,'linear');
% mesh(xq*L_r/2,yq*L_q/2,LvN_f(xq*L_r/2, yq*L_q/2));
plot3(xq*L_r/2,yq*L_q/2,LvN_f(xq*L_r/2, yq*L_q/2), '.');
figure(figurecounter)
figurecounter = figurecounter + 1;
plot3(r_coordinates, q_coordinates,pot_vector, '.');
% mesh(xq,yq,vq);

%% set up System Matrix
[A,M] = LvN_IPDG2D(pot); % Setup using LvN_IPDG2D.m

%% set up Dirichlet boundary conditions
uD = zeros(Nfp*Nfaces, K);
uD(mapD) = exp(-Fx(mapD)).*exp(-Fy(mapD));
% divide mapD into q=+-L_q/2 and r=+-L_r/2 parts
mapSides =  find(abs(abs(Fx)-max(Fx(:))) < 1e-12);
mapTopBot = find(abs(abs(Fy)-max(Fy(:))) < 1e-12 & abs(abs(Fx)-max(Fx(:))) > 1e-12);
% extract only those nodes, which belong to mapD (this removes nodes which
% are doubled as they belong to two different simplices (the corners)
%drin_s = ismember(mapSides,mapD);
%drin_tb = ismember(mapTopBot,mapD);
%mapTopBot = mapTopBot(drin_tb(:));
%mapSides = mapSides(drin_s(:));

% Calculate mu
%q_f = Fy(mapSides) * L_q / 2; % 1d real coordinates at the two sides left and right
q_f = Fy(mapSides);
%q_f = q_f(find(q_f>0));
mu = newtonRaphson(@nullstellenSucheMu, 1.5*e);



fermDiracFt = @(k_value) fermi_dirac_ft(k_value, q_f, mu);
upper_k = sqrt(20*m*mu/hbar/hbar);
f_hut = 2/(2*pi)*integral(fermDiracFt, 0, upper_k, 'ArrayValued', true);

% f_hut_trapez = zeros(size(q_f,1),1);
% safety_k = 1.3;
% nr_steps_k = 5e6;
% steplength_k = safety_k*upper_k/nr_steps_k;
% 
% k = linspace(0 , safety_k*upper_k , nr_steps_k);
% for i=1:size(q_f,1)
% %     f_hut(i) = 2/(2*pi)*integral( @(k)fermDiracFt(k, q(i)),0,2e10);
% %     f_hut(N_q-i+1) = f_hut(i);
%     y_temp = fermDiracFt(k, q_f(i));
%     f_hut_trapez(i) = steplength_k * 2/(2*pi)*trapz(y_temp);
% end
% clear k y_temp;

uD(mapSides) = f_hut;
uD(mapTopBot) = 0;

% evaluate boundary condition contribution to rhs
qN = zeros(Nfp*Nfaces, K);
Aqbc = LvN_IPDGbc2D(uD, qN);
%%

% set up right hand side forcing
rhs = Aqbc;
%rhs = 2*exp(-x).*exp(-y);
%rhs = -2*(pi^2)*sin(pi*x).*sin(pi*y);
%rhs = -MassMatrix*(J.*rhs) + Aqbc;

%% solve system
u = (A-M)\rhs(:);
%u = A\rhs(:);
u = reshape(u, Np, K);

%[rq,tq] = meshgrid(0:0.01:1, 0:pi/40:2*pi);
%xq = rq.*cos(tq);
%yq = rq.*sin(tq);
%%%%% plot result
figure(figurecounter)
figurecounter = figurecounter + 1;
[xq,yq] = meshgrid(-1:0.01:1, -1:0.02:1);
vq = griddata(x,y,u,xq,yq,'cubic');
%vq = griddata(x,y,u,xq,yq,'cubic');
mesh(xq,yq,vq);
%plot3(x,y,u,'o');