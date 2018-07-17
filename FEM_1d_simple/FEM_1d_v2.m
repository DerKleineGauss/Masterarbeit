%%%% Initialisierungen %%%%
% Variablen
N_q = 5e1;     % Anzahl Diskretisierungspunkte in q-Richtung
N_r = 15e1;     % Anzahl Diskretisierungspunkte in r-Richtung
L_q = 152e-9;       % nm, lukas paper
L_r = 106e-9;       % nm, lukas paper
Delta_q = L_q / N_q;
Delta_r = L_r / N_r;
global m kb e hbar Temperature beta N_D E_c;
m=9.1e-31 *0.063;              % Elektronenmasse
kb=1.38064852e-23;
Temperature=300;
beta = 1/(kb*Temperature);
hbar=6.626070040e-34/(2*pi);
e = 1.60217662e-19;     % e-charge
E_c = 1.424*e;          % Leitungsbandkante von GaAs (Valenzbandkante auf 0 gesetzt)
N_D = 1e22;             % Donatorkonzentration
% Calculate mu
%mu = fzero(@(mu)nullstellenSucheMu(mu), [0,3*e]);
mu = newtonRaphson(@nullstellenSucheMu, 1.5*e);

% Rechengebiet
r = zeros(N_r,1);
q = zeros(N_q,1);
for n=1:N_r
    r(n) = (n-1/2)*Delta_r - L_r/2;
end
for k=1:N_q
    q(k) = (k-1/2)*Delta_q - L_q/2;    % 1 ... N-1
end
r_interface = zeros(N_r+1,1);
q_interface = zeros(N_q+1,1);
for n=1:N_r+1
    if (n==N_r+1)
        r_interface(n) = r(n-1)+Delta_r/2;
    else
        r_interface(n) = r(n)-Delta_r/2;
    end
end
for k=1:N_q+1
    if (k==N_q+1)
        q_interface(k) = q(k-1)+Delta_q/2;
    else
        q_interface(k) = q(k)-Delta_q/2;
    end
end

%%
%%%%%% build A %%%%%
alpha = 1/2;
A = full(gallery('tridiag',(N_q),-1,0,1))*alpha;
dim_A = size(A(1,:),2);
a = 0;
b = -1;
c = 1;
b_times_c = b*c;
factor = sqrt(b_times_c);
% sqrt(bc) = sqrt(b/c) = 1i
eigs_A = zeros(dim_A,1);
Q_invers = zeros(dim_A);
norm_Q=0;
dim_A = size(A(1,:),2);
for k=1:dim_A
    eigs_A(k) = (a - 2*factor*cos(pi*k/(dim_A+1)))*alpha;
    for n=1:dim_A
        Q_invers(k,n) = (1i^n)*sin(n*pi*k/(dim_A+1));
    end
    norm_Q = norm_Q+Q_invers(k,1)^2;
end
Q_invers = 1/sqrt(norm_Q) * Q_invers;
Q = Q_invers';
%eigs_A = flipud(eigs_A);
Lambda = diag(eigs_A);
Lambda_invers = zeros(dim_A);
for k=1:dim_A
    Lambda_invers(k,k)=1/Lambda(k,k);
end
% Tests
%norm(eye(dim_A)-Q'*Q)
%Test = (A + Q'*Lambda_A*Q);
%norm(Test)
%norm((A-eigs_A(1)*eye(dim_A))*Q(:,1))
%%

%%%%%% build L %%%%%
% init quadratic L, the System Matrix of size #firstArgument x #secondArgument 
% with #thirdArgument nonzeros
L = spalloc( N_r*N_q , N_r*N_q , 2*N_r*N_q*N_q);    % size: 2 "traces": T and F with each Nq x Nq matrices
%L = zeros(N_r*N_q);
% prepare T,F
[T,F] = deal(zeros(N_q), zeros(N_q));
[F_plus,F_minus,T_plus,T_minus] = deal(zeros(N_q, N_q/2), zeros(N_q, N_q/2), zeros(N_q, N_q/2), zeros(N_q, N_q/2));


% Randbedingungen
rho_l = 0;  % not further needed here, should be considered in calculating C (get_C)
rho_r = 0;  % not further needed here, should be considered in calculating C (get_C)
f_hut = zeros(size(q,1),1);
f_hut_trapez = zeros(size(q,1),1);
fermDiracFt = @(k_value, q_value) fermi_dirac_ft(k_value, q_value, mu);
k = linspace(0,2e10,5e6);
for i=1:N_q/2
%     f_hut(i) = 2/(2*pi)*integral( @(k)fermDiracFt(k, q(i)),0,2e10);
%     f_hut(N_q-i+1) = f_hut(i);
    y = fermDiracFt(k, q(i));
    f_hut_trapez(i) = 2/(2*pi)*trapz(y,k);
    f_hut_trapez(N_q-i+1) = f_hut_trapez(i);
end
clear k y;

startingWithNegativeEV = false;
if imag(eigs_A(1)) < 0
    startingWithNegativeEV = true;
end
if startingWithNegativeEV
   Psi_p = Q_invers(N_q/2+1:N_q ,:) * f_hut_trapez;
   Psi_n = Q_invers(1:N_q/2 , :) * f_hut_trapez;
else
   Psi_n = Q_invers(N_q/2+1:N_q ,:) * f_hut_trapez;
   Psi_p = Q_invers(1:N_q/2 , :) * f_hut_trapez;  
end
b2 = zeros(N_q*N_r,1);

% r-Schleife
for n=1:N_r
    %%%%%% build \Gamma %%%%%
    Gamma_i  = Q_invers * get_C(r_interface(n), N_q, q_interface, Delta_q) * Q;
    Gamma_i_next = Q_invers * get_C(r_interface(n+1), N_q, q_interface, Delta_q) * Q;

    T = 1i*hbar/(m*Delta_r)*Lambda - 1i*e/(2*hbar)*Gamma_i_next;
    F = 1i*hbar/(m*Delta_r)*Lambda + 1i*e/(2*hbar)*Gamma_i;
    if (n==1 || n==N_r)
        if not(startingWithNegativeEV)
            T_minus =   T(:,1:(N_q/2));
            F_minus =   F(:,1:(N_q/2));
            T_plus =    T(:,(N_q/2 + 1):N_q);
            F_plus =    F(:,(N_q/2 + 1):N_q);
        else
            T_plus =    T(:,1:(N_q/2));
            F_plus =    F(:,1:(N_q/2));
            T_minus =   T(:,(N_q/2 + 1):N_q);
            F_minus =   F(:,(N_q/2 + 1):N_q);
        end
        if (n==1)
            L(1:N_q, 1:N_q/2) = -F_minus;   % Spalten 1 bis Nq/2
            L(1:N_q, N_q/2+1 : N_q*3/2) = T;
            b_2(1:N_q) = -F_plus*Psi_p;
        else
            % Test
            % n*N_q - ((n-1/2)*N_q+1) - N_q/2 % should be -1 ( intervall 1:2 has 2 entries but 2-1=1 ...)
            L( (n-1)*N_q+1 : n*N_q , (n-1/2)*N_q+1 : n*N_q ) = T_plus;    % Spalten n*N_q - N_q/2 +1 bis n*N_q
            L( (n-1)*N_q+1 : n*N_q , (n-3/2)*N_q+1 : (n-1/2)*N_q ) = -F;
            b_2( (n-1)*N_q+1 : n*N_q ) = T_minus*Psi_n;
        end
    else
        L( (n-1)*N_q+1 : n*N_q , (n-3/2)*N_q+1 : (n-1/2)*N_q ) = -F;
        L( (n-1)*N_q+1 : n*N_q , (n-1/2)*N_q+1 : (n+1/2)*N_q ) =  T;
    end
end

%%
    % Solve LGS
Psi = L\transpose(b_2);
    % transform back
% this one works well but uses much memory
% rho2 = kron(eye(N_r), Q_invers)*Psi;
rho = zeros(N_r*N_q,1);
for i=1:N_r
   rho( (i-1)*N_q + 1 : i*N_q ) = Q_invers*Psi( (i-1)*N_q + 1 : i*N_q );
end


%%
%%%% Plotting %%%%

close all

% L
% figure('Name', 'Matrix L -- real part');
% x=1:N_q*N_r;
% y=1:N_q*N_r;
% [X,Y] = meshgrid(x,y);
% mesh(X,Y,abs(real(L)));
% figure('Name', 'Matrix L -- imag part');
% mesh(X,Y,abs(imag(L)));
if N_r*N_q < 1e3
    figure('Name', 'Matrix L -- real part');
    surf(real(L));
    figure('Name', 'Matrix L -- imag part');
    surf(imag(L));
end

% rho
rho = reshape(rho, [N_q,N_r]);
figure('Name', 'rho -- real part');
surf(r,q,real(rho));
view(2);   
figure('Name', 'rho -- imag part');
surf(r,q,imag(rho));
view(2);  


%%


