%%
% pos = permute(pos, [3,2,1]);
% pos_s_2 = reshape(pos, [1 size(pos)]) / 1.6;

pos = permute(pos_tracks, [3, 2, 1]);
pos_c = reshape(pos, [1, size(pos)]);

pos = permute(pos_tracks_s, [3, 2, 1]);
pos_s = reshape(pos, [1, size(pos)]) / (1000/340.5); % Since this is already scaled
pos_s = pos_s(1, 1:78, :, :);

% pos_s = cat(4, pos_s, nan(1, 78, 2, 250));
% pos_s = cat(2, pos_s, pos_s_2);

%% Interpolate the missing data...
% Interpolate single-cell data
% pos = nan(841, 2, 105);
% for i = 1:105
%     x = pos_tracks_s(:, 1, i);
%     y = pos_tracks_s(:, 2, i);
%     ind = ~isnan(x) & ~isnan(y);
%     x = x(ind);
%     y = y(ind);
%     t = 1:sum(ind);
%     fx = fit(t', x, 'smoothingspline', 'smoothingparam', 0.9);
%     fy = fit(t', y, 'smoothingspline', 'smoothingparam', 0.9);
%     t_interp = 1:0.2:sum(ind);
%     x_interp = fx(t_interp);
%     y_interp = fy(t_interp);
%     pos(1:length(t_interp), 1, i) = x_interp;
%     pos(1:length(t_interp), 2, i) = y_interp;
% end
% pos = permute(pos, [3, 2, 1]);
% pos_s = reshape(pos, [1, size(pos)]) / (1000/340.5);
    
%%
c = CellCluster(pos_c);
s = CellCluster(pos_s);

msd_c = Msd(c);
msd_s = Msd(s);

% auto_correl_c = Auto_correl_cme(c);
% auto_correl_s = Auto_correl_cme(s);

auto_correl_c = Auto_correl(c);
auto_correl_s = Auto_correl(s);

pLen_c = cme_len(c);
pLen_s = cme_len(s);

pTime_c = cme_time(c);
pTime_s = cme_time(s);

rect = [50 50 600 500];

% Auto-Correlation
ms = 200; lw = 3; malpha = 0.3;
new_fig([50 50 650 550], 4.5, 45, [0 5], [-0.1 0.5], [], [], "");
plot_scatter(c, ms, lw, malpha, c.blu, (0:c.T-2)/30, nanmean(auto_correl_c, 1), ...
    nanstd(auto_correl_c, 0, 1));
plot_scatter(s, ms, lw, malpha, c.ora, (0:s.T-2)/30, nanmean(auto_correl_s, 1), ...
    nanstd(auto_correl_s, 0, 1));
% eb = shadedErrorBar((0:c.T-2)/30, nanmean(auto_correl_c, 1)', nanstd(auto_correl_c, 0, 1));
% eb.patch.FaceColor = c.blu; eb.edge(1).LineStyle = 'none'; eb.edge(2).LineStyle = 'none';
% eb = shadedErrorBar((0:s.T-2)/30, nanmean(auto_correl_s, 1)', nanstd(auto_correl_s, 0, 1));
% eb.patch.FaceColor = s.ora; eb.edge(1).LineStyle = 'none'; eb.edge(2).LineStyle = 'none';
% plot((0:c.T-2)/30, nanmean(auto_correl_c, 1), 'linewidth', lw, 'color', c.blu);
% plot((0:s.T-2)/30, nanmean(auto_correl_s, 1), 'linewidth', lw, 'color', c.ora);
% scatter((0:c.T-2)/30, nanmean(auto_correl_c, 1), ms, 'markerfacecolor', c.blu, 'markerfacealpha', malpha, 'markeredgealpha', 0);
% scatter((0:s.T-2)/30, nanmean(auto_correl_s, 1), ms, 'markerfacecolor', c.ora, 'markerfacealpha', malpha, 'markeredgealpha', 0);
% Insert fitting line
y = nanmean(auto_correl_c, 1); y = y(1:100);
x = 0:c.T-2; x = x / 30; x = x(1:100);
f1 = fit(x', y', 'exp2');
pl = plot(f1); pl.LineWidth = 4; pl.LineStyle = '-'; pl.Color = [c.blu*0.8 0.8];
y = nanmean(auto_correl_s, 1); y = y(1:100);
x = 0:s.T-2; x = x / 30; x = x(1:100);
f2 = fit(x', y', 'exp2');
pl = plot(f2); pl.LineWidth = 4; pl.LineStyle = '-'; pl.Color = [c.ora*0.8 0.8];
% other settings
r = refline(0, 0); r.LineStyle = ':'; r.LineWidth = 2; r.Color = 'k';
xt = 0:2:4; xticks(xt);
xt = string(xt*60); xticklabels(xt);
yt = 0:0.2:0.4; yticks(yt)
ylabel('\itC (t)')
xlabel('t (min)')
setfigprop(800,800,4, 40)
legend off

% Msd
ms = 100; lw = 3; malpha = 0.3;
yc = msd_c;
ys = msd_s;
xc = 120*(0:c.T-1); % Seconds
xs = 120*(0:s.T-1); % Seconds
% xc = 2*(0:c.T-1); % Minutes
% xs = 2*(0:s.T-1); % Minutes
yc = yc(:, 1:5:end); ys = ys(:, 1:5:end); xc = xc(:, 1:5:end); xs = xs(:, 1:5:end);
new_fig(rect, 4.5, 45, [], [], [], [], "");
plot_scatter(c, ms, lw, malpha, c.blu, log10(xc), log10(nanmean(yc, 1)), ...
    nanstd(log10(yc), 0, 1));
plot_scatter(s, ms, lw, malpha, c.ora, log10(xs), log10(nanmean(ys, 1)), ...
    nanstd(log10(ys), 0, 1));
% eb = shadedErrorBar(log10(xs), log10(nanmean(ys, 1)), ...
%     nanstd(log10(ys), 0, 1)); eb.patch.FaceColor = c.ora;
% eb = shadedErrorBar(log10(xc), log10(nanmean(yc, 1)), ...
%     nanstd(log10(yc), 0, 1)); eb.patch.FaceColor = c.blu;
% plot(log10(xc), log10(nanmean(yc, 1)), 'linewidth', lw, 'color', c.blu); hold on
% plot(log10(xs), log10(nanmean(ys, 1)), 'linewidth', lw, 'color', c.ora);
% scatter(log10(xc), log10(nanmean(yc, 1)), ms , 'markerfacecolor', c.blu,'markeredgealpha',0, 'markerfacealpha', malpha); hold on
% scatter(log10(xs), log10(nanmean(ys, 1)), ms , 'markerfacecolor', c.ora,'markeredgealpha',0, 'markerfacealpha', malpha);
xt = 1:5; xticks(xt);
xt = xt - log10(60); xticklabels(round(xt,1)); % Convert to minute
xlabel("log t (min)")
ylabel("log \delta^2 (\mum^2)")
% xtl = ["10^{3}", "10^4", "10^{5}"]; xticklabels(xtl)
yt = 1:5; yticks(yt)
% ytl = ["10^2", "10^3", "10^4", "10^5"]; yticklabels(ytl)
% xlim([3 4.93])
% xlim([3 5])
% xlim([2.9 5])
% ylim([2 5])
% xlim([1.2218 3.2218])
% r = refline(1, -1.05); r.LineStyle = '--'; r.LineWidth = 1.5; r.Color = c.gra;
% r = refline(2, -3.85); r.LineStyle = '--'; r.LineWidth = 1.5; r.Color = c.gra;
r = line([2.8, 3.5], [1.7, 3.1]); r.LineStyle = '--'; r.LineWidth = 2.5; r.Color = c.gra;
r = line([4.0, 4.7], [3.3, 4.0]); r.LineStyle = '--'; r.LineWidth = 2.5; r.Color = c.gra;
xlim([3 5])
ylim([2 5])
setfigprop(800,800,4, 40)
% fill([log10(xs(t0)) 4.57 4.57 log10(xs(t0))], [3.7 3.7 4.2 4.2], [0.15,0.15,0.15 0.15]);
f = fill([4.0 4.8 4.8 4.0], [3.5 3.5 4.4 4.4], 'r'); f.FaceColor = 'none'; f.LineWidth = 4; f.EdgeAlpha = 0.35;

% Msd inset
ms = 70; lw = 4; malpha = 0.5;
new_fig([102,628,825,514], 4.5, 45, [], [], [], [], "");
plot(log10(xc), log10(nanmean(yc, 1)), 'linewidth', lw, 'color', c.blu); hold on
plot(log10(xs), log10(nanmean(ys, 1)), 'linewidth', lw, 'color', s.ora);
scatter(log10(xc), log10(nanmean(yc, 1)), ms , 'markerfacecolor', c.blu,'markeredgealpha',0, 'markerfacealpha', malpha); hold on
scatter(log10(xs), log10(nanmean(ys, 1)), ms , 'markerfacecolor', s.ora,'markeredgealpha',0, 'markerfacealpha', malpha);
% plot(xc, nanmean(yc, 1), 'linewidth', lw, 'color', c.blu); hold on
% plot(xs, nanmean(ys, 1), 'linewidth', lw, 'color', s.ora);
% scatter(xc, nanmean(yc, 1), ms , 'markerfacecolor', c.blu,'markeredgealpha',0, 'markerfacealpha', malpha); hold on
% scatter(xs, nanmean(ys, 1), ms , 'markerfacecolor', s.ora,'markeredgealpha',0, 'markerfacealpha', malpha);
% set(gca, 'yscale', 'log', 'xscale', 'log')
% xt = 0:2000:84000; xticks(xt);
% xtl = repmat("", size(xt)); xtl(11) = "4"; 
% xticklabels(xtl);
% xlim([10^3.8 10^4.7])
% ylim([10^3.2 10^4.4])
% xlabel("t (sec)")
% ylabel("\delta^2 (\mum^2)")
xt = [4 4.5]; xticks(xt)
% xtl = ["10^4", "10^{4.5}"]; xticklabels(xtl)
% xtl = ["4", "4.5"]; xticklabels(xtl);
yt = 3.2:0.6:4.8; yticks(yt)
% ytl = "10^{" + string(yt) + "}"; yticklabels(ytl)
% hold on; r = refline(alpha1 , log10(myc(t0)) - alpha1 * log10(xc(t0))); r.LineWidth = 2; r.Color = 'k'; r.LineStyle = '-';
% hold on; r = refline(alpha2 , log10(mys(t0)) - alpha2 * log10(xs(t0))); r.LineWidth = 2; r.Color = 'k'; r.LineStyle = '-';
% xlim([log10(xs(t0)) 4.57])
% xlim([3.8 4.8])
% ylim([3.2 4.4])
xlim([2.2+log10(60) 4.8])
ylim([3.5 4.4])
xt = (2.2+log10(60)):0.4:(3.0+log10(60)); xticks(xt)
xt = xt - log10(60); xticklabels(round(xt,3)); % Convert to minute
% xlabel("log_{10} t (min)")
ylabel("log \delta^2 (\mum^2)")
r = line([4.3, 4.7], [3.6, 4.0]); r.LineStyle = '--'; r.LineWidth = 3.5; r.Color = c.gra;
% scatter(2.4 + log10(60), 3.6, 300, '^', 'k', 'filled');
% scatter(2.7 + log10(60), 3.9, 300, '^', 'k', 'filled');
scatter(4.202, 0.8, 300, '^', 'k', 'filled');
scatter(4.525, 0.8, 300, '^', 'k', 'filled');
% ylim([3.7 4.2])
% xlim([3.0 4.7])
setfigprop(600,350,4, 40)

% Alpha (do this first than Msd inset)
lw = 4; ms = 130; malpha = 0.3;
new_fig([102,628,825,514], 4.5, 45, [], [2.7 5.3], [], [], "");
myc = squeeze(nanmean(yc, 1));
mys = squeeze(nanmean(ys, 1));
% dt = 1;
dt = 10;
numc = log10(myc(1+dt:end)) - log10(myc(1:end-dt));
dnumc = log10(xc(1+dt:end)) - log10(xc(1:end-dt));
alphac = numc ./ dnumc;
nums = log10(mys(1+dt:end)) - log10(mys(1:end-dt));
dnums = log10(xs(1+dt:end)) - log10(xs(1:end-dt));
alphas = nums ./ dnums;
% yyaxis left
new_logxc = (log10(xc(1:end-dt)) + log10(xc(1+dt:end)))/2;
new_logxs = (log10(xs(1:end-dt)) + log10(xs(1+dt:end)))/2;
plot(new_logxc, alphac, '-', 'linewidth', lw, 'color', [c.blu]);
plot(new_logxs, alphas, '-', 'linewidth', lw, 'color', [c.ora]);
scatter(new_logxc, alphac, ms , 'markerfacecolor', c.blu,'markeredgealpha',0, 'markerfacealpha', malpha); hold on
scatter(new_logxs, alphas, ms , 'markerfacecolor', s.ora,'markeredgealpha',0, 'markerfacealpha', malpha);
% xlim([2.7 4.8])
% ylim([-1 2])
ylim([0.5 1.5])
% xlim([3.8 4.8])
xlim([4.0 4.8])
% xt = (1.2+log10(60)):0.4:(2.8+log10(60)); xticks(xt)
xt = (2.2+log10(60)):0.4:(3.0+log10(60)); xticks(xt)
xt = xt - log10(60); xticklabels(round(xt,3)); % Convert to minute
r = line([0 5], [1 1]); r.LineWidth = 2; r.LineStyle = '--'; r.Color = 'k';
% xlabel("log_{10} t (min)")
ylabel("\alpha")
% scatter(2.4 + log10(60), 0.8, 300, '^', 'k', 'filled');
% scatter(2.7 + log10(60), 0.8, 300, '^', 'k', 'filled');
scatter(4.202, 0.8, 300, '^', 'k', 'filled');
scatter(4.525, 0.8, 300, '^', 'k', 'filled');
xlim([2.2+log10(60) 4.8])
setfigprop(600,350,4, 40)


yyaxis right
plot(log10(xc), log10(myc), 'x-', 'linewidth', 4, 'color', c.colors(1, :));
plot(log10(xs), log10(mys), 'x-', 'linewidth', 4, 'color', c.colors(2, :));
% Linear slope
t1 = 45; t0 = 26;
alpha1 = (log10(myc(t1)) - log10(myc(t0))) / (log10(xc(t1)) - log10(xc(t0)));
hold on; r = refline(alpha1 , log10(myc(t0)) - alpha1 * log10(xc(t0))); r.LineWidth = 1;
r.Color = 'k';
t1 = 45; t0 = 26;
alpha2 = (log10(mys(t1)) - log10(mys(t0))) / (log10(xs(t1)) - log10(xs(t0)));
hold on; r = refline(alpha2 , log10(mys(t0)) - alpha2 * log10(xs(t0))); r.LineWidth = 1;
r.Color = 'k';

% Shape index 
% Import shape index from C:\Users\HyunGyu\Desktop\Weka segmentation\si(1).csv
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(shapei, 'binwidth', 0.3, 'normalization', 'probability');
h.FaceColor = 0.01*[1 1 1]; h.EdgeAlpha = 0.1; h.LineWidth = 1.5; h.EdgeColor = 'w';
h.FaceAlpha = 0.7;
xlabel('\itp')
% xticks(120:120:480)
% xticks(3:2:11)
xticks(4:2:12)
xlim([3.4 10.5])
scatter(nanmean(shapei), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(630, 630, 4, 40)

% Speed fluctuation
v_c = diff(pos_c, 1, 4);
v_c = absvec(v_c, 3);
v_c = squeeze(v_c) / 2 * c.scale;  % um/min unit
a = v_c([1,2,3,4,5,6], :); 

%% 3D trajectory
t0 = 1;
n0 = 24;
nind = neigh(t0, n0, :);
nind = squeeze(nind);
nind = nind(~isnan(nind));
posn = pos(:, :, nind);
sp = 0.005;
for i=1:length(nind)
    for dim=1:2
        x = posn(:,dim,i);
        xCopy = x(~isnan(x));
        t = 0:length(xCopy)-1; t = t';
        f = fit(t,xCopy,'smoothingspline','SmoothingParam',sp);
        posn(~isnan(x), dim, i) = f(t);
    end
end

rect = [50, 50, 700, 900];

new_fig([50 50 600 900], 2.5, 25, [], [], [], [], ""); hold on 
viewangle = [1496, 1969, 47];
timespan = 1:701;

plot3d_traj(c, nind(1), 1:701, rect, 200, viewangle, c.colors(1,:));
plot3d_traj(c, nind(2), 1:701, rect, 200, viewangle, c.colors(2,:));
plot3d_traj(c, nind(3), 1:701, rect, 200, viewangle, c.colors(3,:));
plot3d_traj(c, nind(4), 1:701, rect, 200, viewangle, c.colors(4,:));
plot3d_traj(c, nind(5), 1:701, rect, 200, viewangle, c.colors(2,:));


%% Save all figures
figHandles = findall(0,'Type','figure'); 

savefig(figHandles, "figures\figs_conf_single_dec_18.fig");






%%
rect = [500 500 700 700];
new_fig(rect, 2.5, 25, [3.5 4.5], [], [], [], "");

% plot(log10(120*(0:c.T-1)), log10(nanmean(msd_c, 1)), 'x'); hold on
% plot(log10(120*(0:s.T-1)), log10(nanmean(msd_s(1:78,:), 1)), 'x');

shadedErrorBar(log10(120*(0:s.T-1)), log10(nanmean(msd_s(1:78,:), 1)), nanstd(log10(msd_s(1:78,:)), 0, 1));
shadedErrorBar(log10(120*(0:c.T-1)), log10(nanmean(msd_c, 1)), nanstd(log10(msd_c), 0, 1));

plot(log10(120*(0:c.T-1)), log10(nanmean(msd_c, 1)), 'x'); hold on
plot(log10(120*(0:s.T-1)), log10(nanmean(msd_s(1:78,:), 1)), 'x');

xtl = {"10^{3.5}", "10^4", "10^{4.5}"};
xticklabels(xtl)
yt = 2:5;
ytl = {"10^2", "10^3", "10^4", "10^5"};
yticks(yt)
yticklabels(ytl)

xlabel("sec")
ylabel("\delta^2 (\mum^2)")

refline(1, -1)
refline(2, -4)


%% Draw alpha with different intervals
ms = 150; lw = 5;
malpha = 0.6;
colors = parula(20);
myc = squeeze(nanmean(yc, 1));
mys = squeeze(nanmean(ys, 1));
new_fig(rect, 4.5, 45, [], [2.7 5.3], [], [], "");
for dt = 10
numc = log10(myc(1+dt:end)) - log10(myc(1:end-dt));
dnumc = log10(xc(1+dt:end)) - log10(xc(1:end-dt));
alphac = numc ./ dnumc;
nums = log10(mys(1+dt:end)) - log10(mys(1:end-dt));
dnums = log10(xs(1+dt:end)) - log10(xs(1:end-dt));
alphas = nums ./ dnums;
% Confluent
new_logxc = (log10(xc(1:end-dt)) + log10(xc(1+dt:end)))/2;
plot(new_logxc, alphac, '-', 'linewidth', lw, 'color', c.blu);
scatter(new_logxc, alphac, ms, 'markerfacecolor', c.blu, 'markerfacealpha', malpha, 'markeredgealpha', 0);
% Single
new_logxs = (log10(xs(1:end-dt)) + log10(xs(1+dt:end)))/2;
plot(new_logxs, alphas, '-', 'linewidth', lw, 'color', c.ora);
scatter(new_logxs, alphas, ms, 'markerfacecolor', c.ora, 'markerfacealpha', malpha, 'markeredgealpha', 0);

% xlim([2.7 4.8])
% xlim([3.4 5])
xlim([3.8 4.7])
ylim([0.5 1.5])
end
r = refline(0, 1); r.Color = 'k';

%% Correlation map
pos = permute(pos_tracks, [3, 2, 1]);
pos = reshape(pos, [1, size(pos)]);
c = CellCluster(pos);

edges = -41:2:41; timespan = 1:c.T-1;
[Mccorrel, Ncount] = correl_map(c, edges, timespan);

figure; 
heatmap(Mccorrel');
% 
% figure;
% heatmap(Ncount');


%% Diamond heatmap (scatter plot)
margin = 16; 
hs = floor((length(edges)+1)/2);

% Mean correlation
M = Mccorrel;
M(hs,hs) = 0;
M = M(hs-margin:hs+margin,hs-margin:hs+margin)';
colorlim = 0.4;
M = M / colorlim;
M(M > 1) = 1;
M(M < -1) = -1;

% Count
% M = Ncount;
% M(hs, hs) = 0;
% M = M(hs-margin:hs+margin,hs-margin:hs+margin)';
% maxcount = 1000;
% M = M / maxcount;

% Draw
figure('pos', [1170,396,517,615]); hold on 
edgesA = edges(1:end-1) + (edges(2) - edges(1))/2;
edgesA = edgesA * 1000 / 340.5;
edgesA = edgesA(hs-margin:hs+margin);
[x, y] = meshgrid(edgesA, edgesA);
x_A = x(1:end-1, 1:end-1) + (edgesA(2) - edgesA(1))/2;
y_A = y(1:end-1, 1:end-1) + (edgesA(2) - edgesA(1))/2;
% figure('pos', [100 100 401 409]); hold on 
% scatter(x(:), y(:), 100 , M(:), 'filled', 'd')
A = (M(1:end-1, 1:end-1) + M(1:end-1, 2:end) + M(2:end, 1:end-1) + M(2:end, 2:end)) / 4;
scatter(x_A(:), y_A(:), 100 , A(:), 'filled', 'd')
scatter(x(:), y(:), 100 , M(:), 'filled', 'd')
colormap(parula)
% colormap(spring)
xt = -80 : 20 : 80;
yt = -80 : 20 : 80;
xticks(xt)
yticks(yt)
set(gca, 'fontsize', 20, 'linewidth', 2)
box on
xlim([-80 80])
ylim([-80 80])
xlabel('x(\mum)')
ylabel('y(\mum)')
set(gcf, 'pos', [1201,414.33,498,593.33])

%% Color bar
figure('pos', [100 100 800 800]); 
cb = colorbar;
cb.FontSize = 55;
cb.TickLength = 0.02;
cb.Ticks = 0:0.25:1;
% cb.TickLabels = string(-0.2:0.1:0.2);
cb.TickLabels = string((0:250:1000)/1000);
cb.LineWidth = 3;
cb.Position = [0.5,0.1,0.05,0.82];
colormap(spring)
axis off

%% Axis
% img1 = imread("img1.png"); % Read saved image (edited in PowerPoint)
% img2 = imread("img2.png");




H = length(img1);

fs = 90; lw = 6;
figure;
imshow(img1);
axis on
set(gca, 'fontsize', 20, 'linewidth', 2)
xt = 1:(H-1)/4:H; xticks(xt);
xtl = string(-80:40:80); xticklabels(xtl);
yt = xt; yticks(yt);
ytl = string(80:-40:-80); yticklabels(ytl);
% xlabel('x(\mum)')
% ylabel('y(\mum)')
set(gca, 'fontsize', fs, 'linewidth', lw)
set(gcf, 'pos', [197,42,1564,1313])
% set(gca, 'outerposition', [0.08,0.13,0.85,0.89]) % To avoid labels being cut when in full screen


figure;
imshow(img2);
axis on
set(gca, 'fontsize', 20, 'linewidth', 2)
xt = 1:(H-1)/4:H; xticks(xt);
xtl = string(-80:40:80); xticklabels(xtl);
yt = xt; yticks(yt);
ytl = string(80:-40:-80); yticklabels(ytl);
% xlabel('x(\mum)')
% ylabel('y(\mum)')
set(gca, 'fontsize', fs, 'linewidth', lw)
set(gcf, 'pos', [197,42,1564,1313])
% set(gca, 'outerposition', [0.08,0.13,0.85,0.89]) % To avoid labels being cut when in full screen

%% Fit auto-correlation to exponential & persistence time, length, velocity (confluent)
% Fit options
lowerA = 0.1; lowerTa = -1;
lowerB = 0.1; lowerTb = -1;
upperA = 1; upperTa = -0.002;
upperB = 1; upperTb = -0.002;

obj = c;
Y = auto_correl_c;

coeffs = [];
index = [];
figure; ii = 1;
for i = 1:obj.N
    x = 0:obj.T-2; x = x(1:150);
    y = Y(i,:); y = y(1:150);
%     x = 0:obj.T-2; x = x(1:100);
%     y = Y(i,:); y = y(1:100);
    ind = find(~isnan(y));
    if length(ind) >= 150
        y = y(ind);
        x = x(ind);
        f = fit(x', y', 'exp2', 'lower', [lowerA lowerTa lowerB lowerTb], 'upper', [upperA upperTa upperB upperTb]);
%         subplot(10, 10, ii)
        ii = ii + 1;
        plot(x,y);
        hold on 
        plot(f); legend off
        hold off
        ylim([-0.2 1])
        coeffs = [coeffs; coeffvalues(f)];
        pause(0.1)
        index = [index ; i];
    end
end

a = - 1 ./ coeffs(:, [2 4]);
a1 = max(a, [], 2);

% figure; histogram(a1, 'binwidth', 10);
a2 = min(a, [], 2);
% figure; histogram(a2);

% Persistence time (exponential decay constant)
plen_c = a1 * 2; % Minutes
outliers = find(plen_c >= 950);
plen_c(outliers) = [];

mplen_c = nanmean(plen_c);
% mplen_s = nanmean(plen_s);
sdplen_c = nanstd(plen_c);
% sdplen_s = nanstd(plen_s);

v_c = diff(pos_c, 1, 4) * c.scale / 2; % um/min
v_c = squeeze(absvec(v_c, 3));
v_c = v_c(index, :);
mv_c = nanmean(v_c, 2);
mv_c(outliers) = [];

% Persistence length (um scale)
plength_c = plen_c .* mv_c; % um unit
mplength_c = nanmean(plength_c);
sdplength_c = nanstd(plength_c);

%% Single 
obj = s;
Y = auto_correl_s;

coeffs = []; index = [];
figure; ii = 1;
for i = 1:obj.N
    x = 0:obj.T-2; x = x(1:150);
    y = Y(i,:); y = y(1:150);
%     x = 0:obj.T-2; x = x(1:100);
%     y = Y(i,:); y = y(1:100);
    ind = find(~isnan(y));
    if length(ind) >= 150
        y = y(ind);
        x = x(ind);
        f = fit(x', y', 'exp2', 'lower', [lowerA lowerTa lowerB lowerTb], 'upper', [upperA upperTa upperB upperTb]);
%         subplot(10, 10, ii)
        ii = ii + 1;
        plot(x,y);
        hold on 
        plot(f); legend off
        hold off
        ylim([-0.2 1])
        coeffs = [coeffs; coeffvalues(f)];
        pause(0.1)
        index = [index ; i];
    end
end

a = - 1 ./ coeffs(:, [2 4]);
a1 = max(a, [], 2);

% figure; histogram(a1, 'binwidth', 10);
% a2 = min(a, [], 2);
% figure; histogram(a2);

plen_s = a1 * 2; % Minutes
outliers = find(plen_s >= 950);
plen_s(outliers) = [];

% mplen_c = nanmean(plen_c);
mplen_s = nanmean(plen_s);
% sdplen_c = nanstd(plen_c);
sdplen_s = nanstd(plen_s);

v_s = diff(pos_s, 1, 4) * s.scale / 2;
v_s = squeeze(absvec(v_s, 3));
v_s = v_s(index, :);
mv_s = nanmean(v_s, 2);
mv_s(outliers) = [];

plength_s = plen_s .* mv_s;
mplength_s = nanmean(plength_s);
sdplength_s = nanstd(plength_s);

%% Check
figure;
subplot(2,1,1)
histogram(plen_c, 'normalization', 'probability', 'binwidth', 30)
subplot(2,1,2)
histogram(plen_s, 'normalization', 'probability',  'binwidth', 30)

%% Velocity
v_c = diff(pos_c, 1, 4) * c.scale / 2; % um/min
v_s = diff(pos_s, 1, 4) * c.scale / 2;
v_c = squeeze(absvec(v_c, 3));
v_s = squeeze(absvec(v_s, 3));
mv_c = nanmean(v_c, 2);
mv_s = nanmean(v_s, 2);

% Conf
49.9647 / 30; % Mean
54.2475 / 30; % Std
% Single
31.1116 / 30; % Mean
46.0614 / 30; % Std

%% Diffusion coefficients (conf)
t = 0:c.T-1;
t = t * 2; % In minites
logt = log10(t);
figure; plot(logt, log10(msd_c)); 

y = log10(msd_c);
dy = y(:, 3) - y(:, 2);
dx = logt(3) - logt(2);
a = dy ./ dx;
y0 = y(:, 2);
x0 = logt(2);
logD = y0 - a .* x0;
D_c = 10.^logD;

%%  Diffusion coefficients (single)
t = 0:s.T-1;
t = t * 2; % In minites
logt = log10(t);
figure; plot(logt, log10(msd_s)); 

y = log10(msd_s);
dy = y(:, 3) - y(:, 2);
dx = logt(3) - logt(2);
a = dy ./ dx;
y0 = y(:, 2);
x0 = logt(2);
logD = y0 - a .* x0;
D_s = 10.^logD;

%% Display histogram
figure('pos', [100 100 500 500]); hold on 
h = histogram(D_s, 'binwidth', 0.25, 'normalization', 'probability'); h.FaceColor = c.ora; h.EdgeColor = 'none';
h = histogram(D_c, 'binwidth', 0.25, 'normalization', 'probability'); h.FaceColor = c.blu; h.EdgeColor = 'none';
xlim([0 10])
setfigprop(630, 630, 4, 40)

