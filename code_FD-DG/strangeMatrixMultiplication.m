function [out] = strangeMatrixMultiplication(R , in)
[inner_dim, outer_dim] = size(in);
out = zeros(inner_dim, outer_dim);
if (inner_dim == 1 || outer_dim == 1)
    error("At least one dimension of input to strangeMatrixMultiplication is 1, but should be greater than 1.");
end
[out_Rx, out_Ry]= size(R);
if (out_Rx ~= out_Ry ||out_Rx ~= outer_dim)
    error("Dimension mismatch in R");
end
for i=1:outer_dim
    for j=1:outer_dim
        out(:,i) = out(:,i) + R(i,j)*in(:,j);
    end
end
out = out(:);
end