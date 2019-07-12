function [G_glob] = LVN_systemmatrix_timeDependant(params, B, R)
    % returns systemmatrix A, right hand side rhs containing boundary
    % conditions and additionally transformation matrix R for diagonalizing
    % the FD scheme. Uses strong formulation.

    % x <-> r
    % y <-> q

    Np = params.Np;
    Ny = params.Ny;       % equals 'N_q' in Lukas Paper : the number of cells
    K = params.K;
    systemsize = K*Np*Ny;
    hy = params.hy;

    % build G (from C (from B))
    G = zeros(K*Np, Ny, Ny);
    for r_index=1:K*Np
        C  =  diag(  B(r_index,1:end-1) + B(r_index,2:end)  ) ...
            + diag(B(r_index,2:end-1), -1) ...
            + diag(B(r_index,2:end-1), +1);
        C = C*hy/4;
        G(r_index,:,:) = R'*C*R * 1i;
    end
    

    %% build global matrix G_gl
    
    % G_glob
    indices_row = zeros(K*Np*Ny*Ny,1); %one diagonal for each y block expected
    indices_col = indices_row;
    values_G_glob = indices_row;
    entries = (1:K*Np)';
    for i = 1:Ny
        rows = (i - 1)*(Np*K) + 1 : i*(Np*K);
        for j = 1:Ny
            cols = (j - 1)*(Np*K) + 1 : j*(Np*K);
            values_G_glob(entries) =  G(:,i,j);
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
    
end
