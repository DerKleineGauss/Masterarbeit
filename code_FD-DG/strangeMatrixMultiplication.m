function [out] = strangeMatrixMultiplication(R , in)
[inner_dim, outer_dim] = size(in);
out = zeros(inner_dim, outer_dim);
for i=1:outer_dim
    for j=1:outer_dim
        out(:,i) = out(:,i) + R(i,j)*in(:,j);
    end
end
end