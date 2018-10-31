function [l] = lagrangePolynoms(x, y, x_node, y_node, i)
l = 1;
for j=1:length(x_node)
    if i~=j
        l = l .* (x-x_node(j)) .* (y-y_node(j)) %./ ( (x_node(i) - x_node(j)) .* (y_node(i)-y_node(j)));
    end
end
