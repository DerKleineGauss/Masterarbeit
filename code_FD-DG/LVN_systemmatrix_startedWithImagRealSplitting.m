function [A] = LVN_systemmatrix(params, B)

    % x <-> r
    % y <-> q
    
    Np = params.Np;
    Npy = params.Npy;   % number of nodes
    Ny = Npy - 1;       % equals 'N_q' in Lukas Paper : the number of cells
    size_y = 2*Ny;      % imag and real part
    K = params.K;
    Nfaces = params.Nfaces;
    Nfp = params.Nfp;
    Fmask = params.Fmask;
    EToE = params.EToE;
    EToF = params.EToF;
    x = params.x;
    NODETOL = params.NODETOL;
    rx = params.rx;
    Dr = params.Dr;
    MassMatrix = params.MassMatrix;
    loc2glb = params.loc2glb;
    J = params.J;
    hy = params.hy;
    
    % matrix caring about time derivation term
    T = zeros(size_y, size_y);
    T(1:Ny , Ny+1:2*Ny) = -hy * eye(Ny);
    T(Ny+1:2*Ny , 1:Ny) = +hy * eye(Ny);
    
    % matrix L standing in front of space derivative
    L = zeros(size_y, size_y);
    a = 0;
    b = -1/2;c = 1/2;b_times_c = b*c;factor = sqrt(b_times_c);
    % sqrt(bc) = sqrt(b/c) = 1i
    eigs_A = zeros(Npy-1,1);
    Q = zeros(Npy-1);
    norm_Q=0;
    for k=1:Npy-1
        eigs_A(k) = (a - 2*factor*cos(pi*k/((Npy-1)+1)));
        for n=1:Npy-1
            Q(k,n) = (b/c)^(n/2)*sin(n*pi*k/Npy);
        end
        norm_Q = norm_Q+Q(k,1)^2;
    end
    Q = 1/sqrt(norm_Q) * Q;

    Lambda = diag(eigs_A);
    Lambda_invers = zeros(Npy-1);
    for k=1:Npy-1
        Lambda_invers(k,k)=1/Lambda(k,k);
    end
    L(1:Ny , Ny+1 : 2*Ny) = +Lambda/1i;
    L(Ny+1 : 2*Ny , 1:Ny) = -Lambda/1i;

        % Tests
        if (params.testing)
            A = full(gallery('tridiag',(Npy-1),-1,0,1)) / 2;
            display(['|| 1 - Q_invers Q || = ', num2str(norm(eye(Npy-1)-Q'*Q))]);
            Test = (A - Q'*Lambda*Q);
            display(['|| A - Q_invers \Lambda Q || = ' , num2str(norm(Test))]);
%             Q_inverse = Q';
%             norm((A-eigs_A(1)*eye(Npy-1))*Q_inverse(:,1))
        end
    
    % next build G (from C (from B))
    G = zeros(K*Np, Ny, Ny); 
    for r_index=1:K*Np
        C  =  diag(  B(r_index,1:end-1) + B(r_index,2:end)  ) ...
            + diag(B(r_index,2:end-1), -1) ...
            + diag(B(r_index,2:end-1), +1);
        G(r_index,:,:) = Q'*C*Q *hy/4;
        G(abs(G)<params.NODETOL) = 0;
    end
    spy(G);
    
    for r_index=1:K
        if(~mod(r_index,1000)) r_index, end;
        rows1 = loc2glb(:,r_index)*ones(1,Np);
        cols1 = rows1';
        % Build local operators
        Dx = rx(1, r_index)*Dr;
        
        OP11 = J(1,r_index)*(Dx + MassMatrix*B(rows1(:,1), cols1(1,:)));
end