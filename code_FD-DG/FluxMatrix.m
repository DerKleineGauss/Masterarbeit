function [flux]= FluxMatrix (params, A, S, Lambda)

    Np = params.Np;
    Npy = params.Npy;   % number of nodes
    Ny = Npy - 1;       % equals 'N_q' in Lukas Paper : the number of cells
    K = params.K;
    systemsize = K*Np*Ny;
    
    loc2glb = 1 : systemsize;
    loc2glb = reshape(loc2glb, Np, K, Ny);
    
    % Definitions from paper sheet 2
    e_l  = params.LIFT(:,1)*params.Fscale(1);   % HERE WE ASSUME AN EQUIDISTANT GRID !!!
    e_r = params.LIFT(:,2)*params.Fscale(1);
    
    absA = S * abs(Lambda) * S';
    a = 0.5 * (A + absA);
    b = 0.5 * (A - absA);
    
    indices_row = zeros(4*Np * K * Ny ,1);
    indices_col = indices_row;
    values = indices_row;
    entries = (1:Np)';
    
    for j=1:Ny  % iterate over rows
        rows = (j - 1)*(Np*K) + 1 : j*(Np*K);
        for m=1:Ny  % iterate over cols
            cols = (m - 1)*(Np*K) + 1 : m*(Np*K);
            for k1=2:K-1
                id_right_minus = loc2glb(Np, k1, m);    % minus : interior information
                id_left_minus  = loc2glb(1, k1, m);     % minus : interior information
                id_left_plus   = loc2glb(Np, k1-1, m);  % plus : exterior information
                id_right_plus  = loc2glb(1, k1+1, m);   % plus : exterior information
                % left minus
                indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
                indices_col(entries) = repmat(id_left_minus, Np, 1);
                values(entries) = - a(j,m) * e_l ;
                entries = entries + Np;
                % left plus
                indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
                indices_col(entries) = repmat(id_left_plus, Np, 1);
                values(entries) =  b(j,m) * e_r ;
                entries = entries + Np;
                % right minus
                indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
                indices_col(entries) = repmat(id_right_minus, Np, 1);
                values(entries) = + a(j,m) * e_r;
                entries = entries + Np;
                % right plus
                indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
                indices_col(entries) = repmat(id_right_plus, Np, 1);
                values(entries) = - b(j,m) * e_l;
                entries = entries + Np;
            end
        end
    end
    flux = [indices_row(:), indices_col(:), values(:)];
    flux = flux(1:max(entries)  -Np,:);
    flux = myspconvert(flux, systemsize, systemsize, 1e-15);