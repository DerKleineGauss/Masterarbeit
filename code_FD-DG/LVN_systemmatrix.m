function [A] = LVN_systemmatrix(params, B)

    % x <-> r
    % y <-> q

    Np = params.Np;
    Npy = params.Npy;   % number of nodes
    Ny = Npy - 1;       % equals 'N_q' in Lukas Paper : the number of cells
    K = params.K;
    systemsize = K*Np*Ny;
    Nfaces = params.Nfaces;
    Nfp = params.Nfp;

    x = params.x;
    NODETOL = params.NODETOL;
    rx = params.rx;
    Dr = params.Dr;
    MassMatrix = params.MassMatrix;
    loc2glb = params.loc2glb;
    J = params.J;
    hy = params.hy;

    % matrix L standing in front of space derivative
    a = 0;
    b = -1/2;c = 1/2;b_times_c = b*c;factor = sqrt(b_times_c);
    % sqrt(bc) = sqrt(b/c) = 1i
    eigs_A = zeros(Npy-1,1);
    Q = zeros(Npy-1);
    norm_Q=0;
    if(mod(Npy,2) == 0)
        error("Number of nodes (cells) must in y direction be odd (even).");
    end
    for k=1:Npy-1
        eigs_A(k) = (a - 2*factor*cos(pi*k/((Npy-1)+1)));
        for n=1:Npy-1
            Q(k,n) = (b/c)^(n/2)*sin(n*pi*k/Npy);
        end
        norm_Q = norm_Q+Q(k,1)^2;
    end
    Q = 1/sqrt(norm_Q) * Q;

    L = diag(eigs_A);
    L_invers = zeros(Npy-1);
    for k=1:Npy-1
        L_invers(k,k)=1/L(k,k);
    end
        % Tests
        if (params.testing)
            A = full(gallery('tridiag',(Npy-1),-1,0,1)) / 2;
            display(['|| 1 - Q_invers Q || = ', num2str(norm(eye(Npy-1)-Q'*Q))]);
            Test = (A - Q'*L*Q);
            display(['|| A - Q_invers \Lambda Q || = ' , num2str(norm(Test))]);
%             Q_inverse = Q';
%             norm((A-eigs_A(1)*eye(Npy-1))*Q_inverse(:,1))
        end

%     [params] = adaptMaps(params, eigs_A);

    % interior face variational terms
    loc2glb = 1 : systemsize;
    loc2glb = reshape(loc2glb, Np, K, Ny);
    lift_left  = params.LIFT(:,1)*params.Fscale(1);
    lift_right = params.LIFT(:,2)*params.Fscale(1);
    for j=1:Ny/2    % lambda has negative imag part
        rows = (j - 1)*(Np*K) + 1 : j*(Np*K);
        cols = rows';
        for k1=2:K-1
            id_right_minus = loc2glb(Np, k1, j);    % minus : interior information
            id_left_minus  = loc2glb(1, k1, j);     % minus : interior information
            id_right_plus  = loc2glb(Np, k1-1, j);  % plus : exterior information
            id_left_plus   = loc2glb(1, k1+1, j);   % plus : exterior information
            % left minus
            indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
            indices_col(entries) = repmat(id_left_minus, Np, 1);
            values(entries) = - eigs_A(j) * lift_left * (0.5 + (1-alpha)/2 * nx(1,1));
            entries = entries + Np;
            % left plus
            indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
            indices_col(entries) = repmat(id_left_plus, Np, 1);
            values(entries) = - eigs_A(j) * lift_left * (0.5 + (1-alpha)/2 * nx(2,1));
            entries = entries + Np;
            % right minus
            indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
            indices_col(entries) = repmat(id_right_minus, Np, 1);
            values(entries) = + eigs_A(j) * lift_right * (0.5 + (1-alpha)/2 * nx(1,1));
            entries = entries + Np;
            % right plus
            indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
            indices_col(entries) = repmat(id_right_plus, Np, 1);
            values(entries) = + eigs_A(j) * lift_right * (0.5 + (1-alpha)/2 * nx(2,1));
            entries = entries + Np;
        end
    end
    for j=Ny/2+1 : Ny   % lambda has positive imag part
        for k1=2:K-1

        end
    end

    % next build G (from C (from B))
    G = zeros(K*Np, Ny, Ny);
    for r_index=1:K*Np
        C  =  diag(  B(r_index,1:end-1) + B(r_index,2:end)  ) ...
            + diag(B(r_index,2:end-1), -1) ...
            + diag(B(r_index,2:end-1), +1);
        G(r_index,:,:) = Q'*C*Q *hy/4;
        G(abs(G)<NODETOL) = 0;
    end

    %% build global matrices T_gl, S_gl, G_gl
    % L
    Dx_all_K = zeros(K*Np, Np);
    for k1=1:K
        Dx_all_K((k1-1)*Np + 1 : k1*Np, :) = rx(1, k1)*Dr;
    end
    % global node numbering
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
            values_L_glob(entries) = L(cell_y,cell_y) * Dx;   % L is Blockdiagonal
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
            values_G_glob(entries) = - G(:,i,j);
            indices_row(entries) = rows;
            indices_col(entries) = cols;
            entries = entries + Np*K;
        end
    end
    G_glob = [indices_row(:), indices_col(:), values_G_glob(:)];
    G_glob = G_glob(1:max(entries)  -Np*K,:);
    G_glob = myspconvert(G_glob, systemsize, systemsize, 1e-15);


    % matrix caring about time derivation term
    T_glob = 1i * hy * speye(systemsize);
    a=1
end
