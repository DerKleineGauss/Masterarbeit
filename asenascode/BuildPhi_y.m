function [Phi, p_DFT] = BuildPhi_y(Kx, Ky, Npx, Npy)

% function [Phi] = BuildPhi(Kx, Ky, Npx, Npy)
% Purpose: Build global Discrete Fourier Transformation matrix with respect
% to y variable

Np = Npx*Npy;
K = Kx*Ky;

M = Ky*Npy;
N_tilde = (M-1)/2;
if (mod(M,2)==1)
    warning('The product M=Kx*Npx is odd, so the p=0 value appears within the DFT');
end
indices = linspace(-N_tilde, N_tilde, M);
p_linear = indices * 2*pi/M;

% p_DFT = zeros(Np,K);
% for k1=1:K
%     if (k1 == Kx)
%         index_low = (Kx-1)*Npx + 1;
%     else
%         index_low = mod(k1-1,Kx)*Npx + 1;
%     end
%     index_high = index_low + Npx - 1;
%   p_DFT(:,k1) = repmat(p_linear(index_low : index_high), 1, Npy); 
% end


% build OP_diag, which will be replicated throughout Phi and which
% represents one row of elements
[n_mesh, p_mesh] = meshgrid(indices, p_linear);
all_exponentials = exp(-1i*p_mesh.*n_mesh);

OP_diag = zeros(M*Npy, K*Np);
OP_diag_nonzeros = zeros(1, K*Np);

counter = 1;
for ky=1:Ky
    for row_x=1:Npx
        OP_diag_nonzeros(: , counter + (row_x-1)*Npy ) = 1;
    end
    counter = counter + Kx*Np;
end
indices = find(OP_diag_nonzeros);
for n=1:Npy
    OP_diag((n-1)*M + 1 : n*M , indices) = all_exponentials;
    % shift
    indices = indices + 1;
end

% assemble Phi
map_nonzeros = find(OP_diag);
constant_row_indices = mod(map_nonzeros, Ky*Np);
constant_row_indices(constant_row_indices==0) = Ky*Np;  % fixing mod problem when number == denominator
constant_col_indices = fix((map_nonzeros-0.1)/(Ky*Np)) + 1; % Ganzzahldivision (aufgerundet, daher die -0.1)

% global node numbering
entries = (1:length(map_nonzeros))';
indices_row = zeros(Kx * length(map_nonzeros),1);
indices_col = zeros(Kx * length(map_nonzeros),1);
values = zeros(Kx * length(map_nonzeros),1);

for kx=1:Kx
    start_index_col = (kx-1)*Np;
    start_index_row = (kx-1)*size(OP_diag, 1);
    indices_row(entries(:)) = start_index_row + constant_row_indices;
    indices_col(entries(:)) = start_index_col + constant_col_indices;
    values(entries(:)) = OP_diag(map_nonzeros);
    entries = entries + length(map_nonzeros);
end

Phi = sparse(indices_row, indices_col, values) / sqrt(M);

spy(Phi);

test = normest(Phi*Phi' - speye(K*Np));
if (test > 1e-14)
    warning('The DFT matrix seems to not be orthogonal. Norm of (Q^T Q - I) is %e', test);
end

return
