function setfigprop(width, height, linewidth, fontsize)
pos = get(gcf, 'pos');
set(gcf, 'pos', [pos(1)/2 pos(2)/2 width height])
set(gca, 'linewidth', linewidth, 'fontsize', fontsize)
end