%% Mean pair correlation
rect = [100 100 800 800];
E = 0:-5:-105;
S = 0:0.5:10.5;
figure('pos', rect); h=heatmap(mc);
% figure('pos', rect); h=heatmap(A);
h.ColorData = h.ColorData(end:-1:1, :);
h.GridVisible = 'off';
h.CellLabelColor = 'none';
% h.ColorLimits = [-0.8,0.2];
% h.XDisplayData = S(1:2:end);
% h.YDisplayData = E(1:2:end);
xl = string(S); 
xl(2:2:end) = ''; xl(3:4:end) = '';
yl = string(E(end:-1:1)); 
yl(2:2:end) = ''; yl(3:4:end) = '';
h.XDisplayLabels = string(xl);
h.YDisplayLabels = string(yl); 
h.Colormap = parula;
h.FontSize = 25;
set(gca,'innerposition', [0 0 1 1])


%% Mean Pearson correlation
figure('pos', [100 100 800 800])
h=heatmap(mpearcor);
h.ColorData = h.ColorData(end:-1:1, :);
h.GridVisible = 'off';
h.CellLabelColor = 'none';
h.ColorLimits = [-0.6,-0.1];
% h.XDisplayData = S(1:2:end);
% h.YDisplayData = E(1:2:end);
xl = string(S); 
xl(2:2:end) = ''; xl(3:4:end) = '';
yl = string(E(end:-1:1)); 
yl(2:2:end) = ''; yl(3:4:end) = '';
h.XDisplayLabels = string(xl);
h.YDisplayLabels = string(yl); 
h.Colormap = cool;
h.FontSize = 25;
set(gca,'innerposition', [0 0 1 1])


%% Mean Distance between constituent cells 
rect = [100 100 800 800];
E = 0:-5:-105;
S = 0:0.5:10.5;
figure('pos', rect); h=heatmap(md * 3);
h.ColorData = h.ColorData(end:-1:1, :);
h.GridVisible = 'off';
h.CellLabelColor = 'none';
h.ColorLimits = [15 60];
% h.XDisplayData = S(1:2:end);
% h.YDisplayData = E(1:2:end);
xl = string(S); 
xl(2:2:end) = '';
yl = string(E(end:-1:1)); 
yl(2:2:end) = '';
h.XDisplayLabels = string(xl);
h.YDisplayLabels = string(yl); 
h.Colormap = hot;
% c = colorbar;
h.FontSize = 25;
set(gca,'innerposition', [0 0 1 1])


%%
T = 10000-1;
t = 0:T; t = t';
logt = log10(t);
dx = diff(logt, 1, 1);

y = log10(msd_th);
dy = diff(y, 1, 1);

alpha_th = dy ./ dx;

M = alpha_th(35:50, :, :);
M = squeeze(nanmean(M, 1));
figure('pos', rect); h=heatmap(M);
h.ColorData = h.ColorData(end:-1:1, :);
h.GridVisible = 'off';
h.CellLabelColor = 'none';
h.ColorLimits = [1 2];
% h.XDisplayData = S(1:2:end);
% h.YDisplayData = E(1:2:end);
xl = string(S); 
xl(2:2:end) = '';
yl = string(E(end:-1:1)); 
yl(2:2:end) = '';
h.XDisplayLabels = string(xl);
h.YDisplayLabels = string(yl); 
h.Colormap = jet;
% c = colorbar;
h.FontSize = 25;
set(gca,'innerposition', [0 0 1 1])


