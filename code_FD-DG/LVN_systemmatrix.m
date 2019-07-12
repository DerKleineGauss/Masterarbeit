function [Systemmatrix, rhs, R] = LVN_systemmatrix(params, B)
    % returns systemmatrix A, right hand side rhs containing boundary
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
    [eigs_A, R] = DiagonalDrift(params);
    
    % interior face variational terms
    [flux, D] = FluxMatrix_strong (params, eigs_A);

    % next build G (from C (from B))
    G = zeros(K*Np, Ny, Ny);
    for r_index=1:K*Np
        C  =  diag(  B(r_index,1:end-1) + B(r_index,2:end)  ) ...
            + diag(B(r_index,2:end-1), -1) ...
            + diag(B(r_index,2:end-1), +1);
        C = C*hy/4;
        G(r_index,:,:) = R'*C*R ;
    end

    %% build global matrices S_gl, G_gl
    
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
            values_L_glob(entries) = eigs_A(cell_y) * Dx;   % L is Blockdiagonal
            indices_row(entries) = rows_inner;
            indices_col(entries) = cols_inner;
            entries = entries + Np*Np;
        end
    end
    L_glob = [indices_row(:), indices_col(:), values_L_glob(:)];
    L_glob = L_glob(1:max(entries)  -Np*Np,:);
    L_glob = myspconvert(L_glob, systemsize, systemsize, 1e-15);

    % G_glob
    indices_row = zeros(K*Np*Ny*Ny,1); %one diagonal for each y block expected
    indices_col = indices_row;
    values_G_glob = indices_row;
    entries = (1:K*Np)';
    for i = 1:Ny
        rows = (i - 1)*(Np*K) + 1 : i*(Np*K);
        for j = 1:Ny
            cols = (j - 1)*(Np*K) + 1 : j*(Np*K);
            values_G_glob(entries) = 1i* G(:,i,j);
            indices_row(entries) = rows;
            indices_col(entries) = cols;
            entries = entries + Np*K;
        end
    end
    G_glob = [indices_row(:), indices_col(:), values_G_glob(:)];
    G_glob = G_glob(1:max(entries)  -Np*K,:);
    G_glob = myspconvert(G_glob, systemsize, systemsize, 1e-15);
    
    if (params.testing)    
        figure('name', "Real part of G_glob");
        spyc_grid(real(G_glob),'cool',Np,K*Np)
        figure('name', "Imag part of G_glob");
        spyc_grid(imag(G_glob),'cool',Np,K*Np)
    end
    
    
    % build rhs consisting of inflow boundary conditions
%     [D_l_pos, D_l_neg, D_r_pos, D_r_neg] = SplitDmatrix(params, D, R, eigs_A);
    [D_l, D_r] = SplitDmatrix(params, D);
    [rhs] = flux_rhs(params, D_l, D_r, R, eigs_A);
    
    % set up systemmatrix A
%     Systemmatrix = L_glob + G_glob - flux - D_l - D_r;    % ist quatsch
%     weil wir keine Werte von außerhalb nehmen, die wir nicht kennen.
    Systemmatrix = L_glob + G_glob - flux ;
    
end
