function plot_sorted(x,y)
labels = cellstr( num2str([1:length(x)]') );  %' # labels correspond to their order

plot(x, y, 'rx', 'markersize', 5)
text(x, y, labels, 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right')