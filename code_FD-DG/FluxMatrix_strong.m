function [flux, D]= FluxMatrix_strong (params, eigs_A)

% Purpose: build global flux matrix for the term 
%   (Lambda u - (Lambda u)*) n M^(-1) \vec(\ell)

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

absLambda = abs(eigs_A);
% wipe out numerical errors
%     ids = find(abs(real(absLambda)) < 1e-15 & abs(real(absLambda)) > 0);
%     absLambda(ids) = 1i*imag(absLambda(ids));
%     ids = find(abs(imag(absLambda)) < 1e-15 & abs(imag(absLambda)) > 0);
%     absLambda(ids) =    real(absLambda(ids));

a = 0.5 * (eigs_A + absLambda);
b = 0.5 * (eigs_A - absLambda);

indices_row = zeros((3*Np*2 + 4*Np*(K-2)) * Ny ,1);
indices_col = indices_row;
values = indices_row;
entries = (1:Np)';

dindices_row = zeros((2*1*Np) * Ny ,1);
dindices_col = dindices_row;
dvalues = dindices_row;
dentries = (1:Np)';

for j=1:Ny  % iterate over rows and cols of block diag matrix
    rows = (j - 1)*(Np*K) + 1 : j*(Np*K);
    
    for k1=1:K
        id_right_minus = loc2glb(Np, k1, j);    % minus : interior information
        id_left_minus  = loc2glb(1, k1, j);     % minus : interior information
        % left minus
        indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
        indices_col(entries) = repmat(id_left_minus, Np, 1);
        values(entries) = - a(j) * e_l ;
        entries = entries + Np;
        % left plus
        if (k1>1)
            id_left_plus   = loc2glb(Np, k1-1, j);  % plus : exterior information
            indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
            indices_col(entries) = repmat(id_left_plus, Np, 1);
            values(entries) = + a(j) * e_l ;
            entries = entries + Np;
        end
        % right minus
        indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
        indices_col(entries) = repmat(id_right_minus, Np, 1);
        values(entries) = + b(j) * e_r;
        entries = entries + Np;
        % right plus
        if (k1<K)
            id_right_plus  = loc2glb(1, k1+1, j);   % plus : exterior information
            indices_row(entries) = rows((k1-1)*Np + 1 : k1*Np);
            indices_col(entries) = repmat(id_right_plus, Np, 1);
            values(entries) = - b(j) * e_r;
            entries = entries + Np;
        end
        
        if (k1 == 1)    % left boundary
            % left plus
            dindices_row(dentries) = rows((k1-1)*Np + 1 : k1*Np);
            dindices_col(dentries) = repmat(id_left_minus, Np, 1);
            dvalues(dentries) = + a(j) * e_l ;
            dentries = dentries + Np;
        elseif(k1 == K) % right boundary
            % right plus
            dindices_row(dentries) = rows((k1-1)*Np + 1 : k1*Np);
            dindices_col(dentries) = repmat(id_right_minus, Np, 1);
            dvalues(dentries) = - b(j) * e_r;
            dentries = dentries + Np;
        end
    end
end
%     ids = find(abs(values(:))>tol);
% A = A(ids, :);
flux = [indices_row(:), indices_col(:), values(:)];
flux = flux(1:max(entries)  -Np,:);
flux = myspconvert(flux, systemsize, systemsize, 1e-15);

D = [dindices_row(:), dindices_col(:), dvalues(:)];
D = D(1:max(dentries)  -Np,:);
D = myspconvert(D, systemsize, systemsize, 1e-15);

%     ids = find(abs(real(flux)) < 1e-15 & abs(real(flux)) > 0);
%     flux(ids) = 1i*imag(flux(ids));
%     ids = find(abs(imag(flux)) < 1e-15 & abs(imag(flux)) > 0);
%     flux(ids) = real(flux(ids));

if (params.testing)
    figure('Name', 'Real part of flux matrix');
    spyc_grid(real(flux),'cool', Np, K*Np);
    %         figure('Name', 'Imag part of flux matrix');
    %         spyc_grid(imag(flux),'cool', Np, K*Np);
    figure('Name', 'Imag part of D matrix');
    spyc_grid(imag(D),'cool', Np, K*Np);
    figure('Name', 'Real part of D matrix');
    spyc_grid(real(D),'cool', Np, K*Np);
end
end