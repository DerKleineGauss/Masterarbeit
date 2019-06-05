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
    % wipe out numerical errors
    ids = find(abs(real(absA)) < 1e-15 & abs(real(absA)) > 0);
    absA(ids) = 1i*imag(absA(ids));
    ids = find(abs(imag(absA)) < 1e-15 & abs(imag(absA)) > 0);
    absA(ids) =    real(absA(ids));
    
    a = 0.5 * (A + absA);
    b = 0.5 * (A - absA);
    
    indices_row = zeros((3*Np*2 + 4*Np*(K-2)) * Ny ,1);
    indices_col = indices_row;
    values = indices_row;
    entries = (1:Np)';
    
    for j=1:Ny  % iterate over rows
        rows = (j - 1)*(Np*K) + 1 : j*(Np*K);
        for m=1:Ny  % iterate over cols
            cols = (m - 1)*(Np*K) + 1 : m*(Np*K);
            for k1=1:K
                id_right_minus = loc2glb(Np, k1, m);    % minus : interior information
                id_left_minus  = loc2glb(1, k1, m);     % minus : interior information
                % left minus
                indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
                indices_col(entries) = repmat(id_left_minus, Np, 1);
                values(entries) = - a(j,m) * e_l ;
                entries = entries + Np;
                % left plus
                if (k1>1)
                    id_left_plus   = loc2glb(Np, k1-1, m);  % plus : exterior information
                    indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
                    indices_col(entries) = repmat(id_left_plus, Np, 1);
                    values(entries) =  b(j,m) * e_r ;
                    entries = entries + Np;
                end
                % right minus
                indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
                indices_col(entries) = repmat(id_right_minus, Np, 1);
                values(entries) = + a(j,m) * e_r;
                entries = entries + Np;
                % right plus
                if (k1<K)
                    id_right_plus  = loc2glb(1, k1+1, m);   % plus : exterior information
                    indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
                    indices_col(entries) = repmat(id_right_plus, Np, 1);
                    values(entries) = - b(j,m) * e_l;
                    entries = entries + Np;
                end
            end
        end
    end
%     ids = find(abs(values(:))>tol);
% A = A(ids, :);
    flux = [indices_row(:), indices_col(:), values(:)];
    flux = flux(1:max(entries)  -Np,:);
    flux = myspconvert(flux, systemsize, systemsize, 1e-15);
    
%     ids = find(abs(real(flux)) < 1e-15 & abs(real(flux)) > 0);
%     flux(ids) = 1i*imag(flux(ids));
%     ids = find(abs(imag(flux)) < 1e-15 & abs(imag(flux)) > 0);
%     flux(ids) = real(flux(ids));
    
    if (params.testing)
        figure('Name', 'Real part of flux matrix');
        spyc_grid(real(flux),'cool', Np, K*Np);
        figure('Name', 'Imag part of flux matrix');
        spyc_grid(imag(flux),'cool', Np, K*Np);
    end