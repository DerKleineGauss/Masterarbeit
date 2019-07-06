function [Systemmatrix, rhs, R, maxSpeed] = LVN_systemmatrix_timeIndependent(params)
    % returns time-independent systemmatrix A (not including potential part), 
    % right hand side rhs containing boundary
    % conditions and additionally transformation matrix R for diagonalizing
    % the FD scheme. Uses strong formulation.

    % x <-> r
    % y <-> q

    Np = params.Np;
    Npy = params.Npy;   % number of nodes
    Ny = params.Ny;       % equals 'N_q' in Lukas Paper : the number of cells
    K = params.K;
    systemsize = K*Np*Ny;
    NODETOL = params.NODETOL;
    rx = params.rx;
    Dr = params.Dr;
    hy = params.hy;

    % matrix A, standing in front of space derivative and its Eigenvalue
    % Decomposition A = S \Lambda S'
    a = 0;
    b = +1i/2;c = -1i/2; factor = sqrt(b*c);
    A = full(gallery('tridiag',(Npy-1),b,0,c)) ;
    % sqrt(bc) = sqrt(b/c) = 1i
    eigs_A = zeros(Npy-1,1);
    R = zeros(Npy-1);
    norm_Q=0;
    if(mod(Npy,2) == 0)
        error("Number of nodes (cells) must in y direction be odd (even).");
    end
    for k=1:Npy-1
        eigs_A(k) = (a + sign(imag(c)) * 2*factor*cos(pi*k/((Npy-1)+1)));
        for n=1:Npy-1
            R(k,n) = (b/c)^(n/2)*sin(n*pi*k/Npy);
        end
        norm_Q = norm_Q+R(k,1)^2;
    end
    R = 1/sqrt(norm_Q) * R;
    R = R';     % we choose R so that A = R \Lambda R' , u=Rv with v being the new variable
    Lambda = diag(eigs_A);
    Lambda_invers = zeros(Npy-1);
    for k=1:Npy-1
        Lambda_invers(k,k)=1/Lambda(k,k);
    end
    maxSpeed = max(abs(eigs_A));
        % Tests
        if (params.testing)
            display(['|| id - R_invers R || = ', num2str(norm(eye(Npy-1)-R'*R))]);
            Test = (A - R*Lambda*R');
            display(['|| A - R \Lambda R_invers || = ' , num2str(norm(Test))]);
            display(['|| (A - \lambda_1 * id)*r_1 || = ', num2str(norm((A-eigs_A(1)*eye(Npy-1))*R(:,1)) )])
        end
    
    % interior face variational terms
    [flux, D] = FluxMatrix_strong (params, Lambda);

    %% build global matrix S_gl
    
    % L
    Dx_all_K = zeros(K*Np, Np);
    for k1=1:K
        Dx_all_K((k1-1)*Np + 1 : k1*Np, :) = rx(1, k1)*Dr;
    end
    indices_row = zeros(K*Np*Ny,1);
    indices_col = indices_row;
    values_L_glob = indices_row;
    entries = (1:Np*Np)';
    for cell_y = 1:Ny
        cols = (cell_y - 1)*(Np*K) + 1 : cell_y*(Np*K);
        rows = cols';

        for k1=1:K
            rows_inner = rows((k1-1)*Np + 1 : k1*Np)*ones(1,Np);
            cols_inner = rows_inner';
            Dx = Dx_all_K((k1-1)*Np + 1 : k1*Np, :);    % = rx(1, k1)*Dr (but saves computational time)
            values_L_glob(entries) = Lambda(cell_y,cell_y) * Dx;   % L is Blockdiagonal
            indices_row(entries) = rows_inner;
            indices_col(entries) = cols_inner;
            entries = entries + Np*Np;
        end
    end
    L_glob = [indices_row(:), indices_col(:), values_L_glob(:)];
    L_glob = L_glob(1:max(entries)  -Np*Np,:);
    L_glob = myspconvert(L_glob, systemsize, systemsize, 1e-15);

    
    % build rhs consisting of inflow boundary conditions
%     [D_l_pos, D_l_neg, D_r_pos, D_r_neg] = SplitDmatrix(params, D, R, eigs_A);
    [D_l, D_r] = SplitDmatrix(params, D);
    [rhs] = flux_rhs(params, D_l, D_r, R, eigs_A);
    
    % set up systemmatrix A
%     Systemmatrix = L_glob + G_glob - flux - D_l - D_r;    % ist quatsch
%     weil wir keine Werte von außerhalb nehmen, die wir nicht kennen.
    Systemmatrix = L_glob + - flux ;
    
end
