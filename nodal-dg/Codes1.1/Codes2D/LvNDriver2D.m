% Driver script for solving the 2D Poisson equation
Globals2D;
GlobalsLvN;
GlobalsLvN_data;

% Polynomial order used for approximation
N = 3;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('block2.neu');

% Initialize solver and construct grid and metric
StartUp2D;

% set up boundary conditions
BuildBCMaps2D;

%% set up right hand side
% evaluate potential term
q_coordinates = y * L_q / 2; % 1d real coordinates
r_coordinates = x * L_r / 2; % 1d real coordinates
pot_vector = LvN_f(r_coordinates,q_coordinates);
% build diagonal matrix representing f(r,q) in the f(r,q)*u part
pot = sparse(1:K*Np, 1:K*Np, pot_vector(:), K*Np, K*Np);
% % plot test
[xq,yq] = meshgrid(-1:0.02:1, -1:0.02:1);
% vq = griddata(x,y,pot_vector,xq,yq,'linear');
mesh(xq,yq,LvN_f(xq*L_r/2, yq*L_q/2));
plot3(xq,yq,LvN_f(xq*L_r/2, yq*L_q/2), '.');
% mesh(xq,yq,vq);


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
q_f = Fy(mapSides) * L_q / 2; % 1d real coordinates at the two sides left and right
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

% solve system
u = (A-M)\rhs(:);
u = reshape(u, Np, K);

%[rq,tq] = meshgrid(0:0.01:1, 0:pi/40:2*pi);
%xq = rq.*cos(tq);
%yq = rq.*sin(tq);
[xq,yq] = meshgrid(-1:0.02:1, -1:0.02:1);
vq = griddata(x,y,u,xq,yq,'cubic');
%vq = griddata(x,y,u,xq,yq,'cubic');
mesh(xq,yq,vq);
%plot3(x,y,u,'o');