% Variablen
N_q = 1e1+1;     % Anzahl Diskretisierungspunkte in q-Richtung
N_r = 1e1;       % Anzahl Diskretisierungspunkte in r-Richtung
L_q = 152;       % nm, lukas paper
L_r = 106;       % nm, lukas paper
Delta_q = L_q / N_q;
Delta_r = L_r / N_r;
m = 1;              % Elektronenmasse
hbar = 1;           % hquer = 1 (n.E.)
e = 1;              % e-charge
alpha = hbar^2/(2*Delta_q*m*e);   % Vorfaktor für die Matrix A (LSchulz_Solution...)

r = zeros(N_r-1,1);
q = zeros(N_q-1,1);
for n=1:N_r-1
    r(n) = n*Delta_r;
end
for k=1:N_q-1
    q(k) = k*Delta_q;    % 1 ... N-1
end

%%%%%% build A %%%%%
A = full(gallery('tridiag',(N_q-1),-1,0,1))*alpha;
a = 0;
b = -1;
c = 1;
b_times_c = b*c;
factor = sqrt(b_times_c);
% sqrt(bc) = sqrt(b/c) = 1i
eigs_A = zeros(N_q-1,1);
Q_invers = zeros(N_q-1);
norm_Q=0;
for k=1:N_q-1
    eigs_A(k) = (a - 2*factor*cos(pi*k/((N_q-1)+1)))*alpha;
    for n=1:N_q-1
        Q_invers(k,n) = (1i^n)*sin(n*pi*k/N_q);
    end
    norm_Q = norm_Q+Q_invers(k,1)^2;
end
Q_invers = 1/sqrt(norm_Q) * Q_invers;
Q = Q_invers';
%eigs_A = flipud(eigs_A);
Lambda = diag(eigs_A);
Lambda_invers = zeros(N_q-1);
for k=1:N_q-1
    Lambda_invers(k,k)=1/Lambda(k,k);
end

% Tests
%norm(eye(N-1)-Q'*Q)
%Test = (A + Q'*Lambda_A*Q);
%norm(Test)
%norm((A-eigs_A(1)*eye(N_q-1))*Q(:,1))

%%%%%% build L %%%%%
B = zeros(N_q-1);
% init non-quadratic L, the System Matrix of size #firstArgument x #secondArgument 
% with #thirdArgument nonzeros
L = spalloc( N_r*(N_q - 1) , (N_r - 1) * (N_q - 1) , 2*(N_r-1)*(N_q-1));

% Randbedingungen
rho_l = 0;
rho_r = 0;


% r-Schleife
for n=1:N_r-1
    % q-Schleife
    for k=1:N_q-1
        % Mittelwertbildung DeltaV_k(r_{n} +- Delta r/2)
        % B(k,k) = 0.5*(delta_V (r(n)+Delta_r/2, q(k), Delta_q) + delta_V(r(n) - Delta_r/2, q(k), Delta_q));
        B(k,k) = delta_V (r(n), q(k), Delta_q) ;
    end
    B(1,1) = B(1,1) - alpha * rho_l;
    B(N_q-1,N_q-1) = B(N_q-1, N_q-1) - alpha * rho_r;

    %%%%%% build \Gamma %%%%%
    Gamma  = Q_invers * B * Q;
    temp = Gamma*Delta_r/2;
    gamma = Lambda - temp;
    delta = Lambda + temp;
    %Gamma_strich = Lambda_A_invers * Gamma ;
    n*(N_q-1)
    if (n<N_r-1)
        % "Diagonalelemente" (L ist nicht quadratisch), dh. die gamma's : erster quadratischer Block zB 1 bis N_q-1
        L( (n-1)*(N_q-1)+1 : n*(N_q-1) , (n-1)*(N_q-1)+1 : n*(N_q-1) ) = gamma;
    end
    if (n>1)
        %L((n-1)*(N_q-1)+1 : n*(N_q-1) , (n-2)*(N_q-1)+1 : (n-1)*(N_q-1)) = delta;
    end
end


