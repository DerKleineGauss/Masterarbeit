function [res] = DFTrhsIndices(Kx, Ky, Npx, Npy, Np)

    const_sorted_x_indices = zeros(Ky*Npy,1);   % indices of points along an x=const line
    res = zeros(Kx*Ky*Np,1);
    entry = 1;
    for ky = 1:Ky
       for npy = 1:Npy
           const_sorted_x_indices(entry) = 1 + (npy-1)*Npx + (ky-1)*Kx*Np;
           entry = entry + 1;
       end
    end
    n = length(const_sorted_x_indices);
    counter = 1;
    for kx = 1:Kx       % gather information in x direction, for each new x=const line 
        % we have to shift the const_sorted_x_indices vector by 1 (still within same cell) or by Np (new cell)
       for npx = 1:Npx
           res(counter : counter+n-1) = const_sorted_x_indices + (npx-1) + (kx-1)*Np;
           counter = counter + n;
       end
    end
end