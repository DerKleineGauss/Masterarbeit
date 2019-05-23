function [Phi, p_DFT] = BuildPhi_y(params)

% function [Phi] = BuildPhi(Kx, Ky, Npx, Npy)
% Purpose: Build global Discrete Fourier Transformation matrix with respect
% to y variable

Np = params.Np;
Npx = params.Npx;
Npy = params.Npy;
K = params.K;
Kx = params. Kx;
Ky = params. Ky;

M = Ky*Npy;
N_tilde = (M-1)/2;
if (mod(M,2)==1)
    warning('The product M=Kx*Npx is odd, so the p=0 value appears within the DFT');
end
indices = linspace(-N_tilde, N_tilde, M);
p_linear = indices * 2*pi/M;

p_DFT = zeros(Np,K);
temp = repmat(p_linear, 1, Npx*Kx);
p_DFT(:) = temp(:);

% prepare row interchanging so that on the rhs there will be u^ (x_1,p?) ...
% u^ (x_last, p_?)
ordering = DFTrhsIndices(Kx, Ky, Npx, Npy, Np);

% check: all points should reside on a vertical line
% figure(82)
% for i=1:Npx*Kx
%     plot_sorted(x(ordering((i-1)*M + 1 : i*M)), y(ordering((i-1)*M + 1 : i*M)));
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

% rearrange (Zeilenvertauschung, damit rechts udach(x_1) bis udach(x_last) sortiert steht
indices_p_DFT = 1:length(p_DFT(:));
indices_p_DFT = ordering(indices_p_DFT);
%p_DFT = p_DFT(indices_p_DFT);
[wayne, I] = sort(indices_p_DFT);
p_DFT = p_DFT(I);
Phi_unsorted = sparse(indices_row, indices_col, values) / sqrt(M);
indices_row = ordering(indices_row);

Phi = sparse(indices_row, indices_col, values) / sqrt(M);

if (params.testing)
    if (sum(abs(Phi_unsorted(2,:) - Phi(1+Npx,:))) > 1e-14)
        warning('AHHHHHHHH');
    end

    % spy(Phi);

    test = normest(Phi*Phi' - speye(K*Np));
    if (test > 1e-14)
        warning('The DFT matrix seems to not be orthogonal. Norm of (Q^T Q - I) is %e', test);
    end
end

return
