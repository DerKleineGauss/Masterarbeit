function [Phi, p_DFT] = BuildPhi_x(Kx, Ky, Npx, Npy)

% function [Phi] = BuildPhi(Kx, Ky, Npx, Npy)
% Purpose: Build global Discrete Fourier Transformation matrix with respect
% to x variable

Np = Npx*Npy;
K = Kx*Ky;

M = Kx*Npx;
N_tilde = (M-1)/2;
if (mod(M,2)==1)
    warning('The product M=Kx*Npx is odd, so the p=0 value appears within the DFT');
end
indices = linspace(-N_tilde, N_tilde, M);
p_linear = indices * 2*pi/M;

p_DFT = zeros(Np,K);
for k1=1:K
    if (k1 == Kx)
        index_low = (Kx-1)*Npx + 1;
    else
        index_low = mod(k1-1,Kx)*Npx + 1;
    end
    index_high = index_low + Npx - 1;
  p_DFT(:,k1) = repmat(p_linear(index_low : index_high), 1, Npy); 
end


% build OP_diag, which will be replicated throughout Phi and which
% represents one row of elements
[n_mesh, p_mesh] = meshgrid(indices, p_linear);
all_exponentials = exp(-1i*p_mesh.*n_mesh);

OP_diag = zeros(M*Npy, Kx*Np);
for kx=1:Kx
    for row_y=1:Npy
        x_l = (kx-1)*Np + (row_y-1)*Npx + 1;
        x_r = x_l + Npx - 1;
        y_l = (row_y-1)*M+1 ;
        y_r =  row_y*M;
        OP_diag(y_l:y_r, x_l : x_r) = all_exponentials(: , (kx-1)*Npx + 1: kx*Npx); % Resultat sieht sinnvoll aus
    end
end

% assemble Phi having Ky times OP_diag on its diagonal
map_nonzeros = find(OP_diag);
constant_row_indices = mod(map_nonzeros, Kx*Np);
constant_row_indices(constant_row_indices==0) = Kx*Np;  % fixing mod problem when number == denominator
constant_col_indices = fix((map_nonzeros-0.1)/(Kx*Np)) + 1; % Ganzzahldivision (aufgerundet, daher die -0.1)

% global node numbering
entries = (1:length(map_nonzeros))';
indices_row = zeros(Ky * length(map_nonzeros),1);
indices_col = zeros(Ky * length(map_nonzeros),1);
values = zeros(Ky * length(map_nonzeros),1);

for ky=1:Ky
    start_index_col = (ky-1)*Np*Kx;
    start_index_row = start_index_col;  % top left element always the start point
    indices_row(entries(:)) = start_index_row + constant_row_indices;
    indices_col(entries(:)) = start_index_col + constant_col_indices;
    values(entries(:)) = OP_diag(map_nonzeros);
    entries = entries + length(map_nonzeros);
end

Phi = sparse(indices_row, indices_col, values) / sqrt(M);

test = normest(Phi*Phi' - speye(K*Np));
if (test > 1e-14)
    warning('The DFT matrix seems to not be orthogonal. Norm of (Q^T Q - I) is %e', test);
end

return
