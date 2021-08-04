%%
% Confluent
i = 1:200;
pos = pos(:, :, i);
pos = permute(pos, [3,2,1]);
pos = reshape(pos, [1 size(pos)]);
cs = CellCluster(pos);
cs.scale = 3; 

% Single cells
i = 1:200;
pos = pos(:, :, i);
pos = permute(pos, [3,2,1]);
pos = reshape(pos, [1 size(pos)]);
ss = CellCluster(pos);
ss.scale = 3;

%%
msd_cs = Msd(cs);
msd_ss = Msd(ss);

% auto_correl_c = Auto_correl_cme(c);
% auto_correl_s = Auto_correl_cme(s);

auto_correl_cs = Auto_correl(cs);
auto_correl_ss = Auto_correl(ss);

pLen_c = cme_len(cs);
pLen_s = cme_len(ss);

pTime_c = cme_time(cs);
pTime_s = cme_time(ss);

rect = [50 50 600 500];

% Auto-Correlation
ms = 150; lw = 6; malpha = 0.8;
new_fig([50 50 650 550], 4.5, 45, [], [], [], [], ""); hold on 
eb = shadedErrorBar((0:cs.T-2)/30 * 60, nanmean(auto_correl_cs, 1)', nanstd(auto_correl_cs, 0, 1)); eb.mainLine.LineStyle = 'none';
eb.patch.FaceColor = cs.blu; eb.edge(1).LineStyle = 'none'; eb.edge(2).LineStyle = 'none';
eb = shadedErrorBar((0:ss.T-2)/30 * 60, nanmean(auto_correl_ss, 1)', nanstd(auto_correl_ss, 0, 1)); eb.mainLine.LineStyle = 'none';
eb.patch.FaceColor = ss.ora; eb.edge(1).LineStyle = 'none'; eb.edge(2).LineStyle = 'none';
% plot_scatter(cs, ms, lw, malpha, cs.blu, (0:cs.T-2)/30, nanmean(auto_correl_cs, 1), ...
%     nanstd(auto_correl_cs, 0, 1));
% plot_scatter(ss, ms, lw, malpha, cs.ora, (0:ss.T-2)/30, nanmean(auto_correl_ss, 1), ...
%     nanstd(auto_correl_ss, 0, 1));
% plot((0:cs.T-2)/30, nanmean(auto_correl_cs, 1), 'linewidth', lw, 'color', cs.blu);
% plot((0:ss.T-2)/30, nanmean(auto_correl_ss, 1), 'linewidth', lw, 'color', cs.ora);
scatter((0:cs.T-2)/30 * 60, nanmean(auto_correl_cs, 1), ms, 'markerfacecolor', cs.blu, 'markerfacealpha', malpha, 'markeredgealpha', 0);
scatter((0:ss.T-2)/30 * 60, nanmean(auto_correl_ss, 1), ms, 'markerfacecolor', cs.ora, 'markerfacealpha', malpha, 'markeredgealpha', 0);
r = refline(0, 0); r.LineStyle = ':'; r.LineWidth = 3; r.Color = 'k';
xt = 0:120:360; xticks(xt)
yt = 0:0.2:0.6; yticks(yt)
% Exponential fits
y = nanmean(auto_correl_cs, 1); y = y(1:180);
x = 0:180-1; x = x * 2;
f = fit(x', y', 'exp2');
pl = plot(x, f(x)); pl.LineWidth = 2; pl.LineStyle = '-'; pl.Color = 'k';
y = nanmean(auto_correl_ss, 1); y = y(1:180);
f = fit(x', y', 'exp2');
pl = plot(x, f(x)); pl.LineWidth = 2; pl.LineStyle = '-'; pl.Color = 'k';
ylabel('\itC (t)')
xlabel('t (min)')
legend off
xlim([0 360])
ylim([-0.05 0.62])
setfigprop(800,800,4, 40)

% Msd
rect = [50 50 600 500];
ms = 80; lw = 3; malpha = 0.8;
yc = msd_cs;
ys = msd_ss;
xc = 120*(0:cs.T-1) / 60;
xs = 120*(0:ss.T-1) / 60;
yc = yc(:, 1:5:end); ys = ys(:, 1:5:end); xc = xc(:, 1:5:end); xs = xs(:, 1:5:end);
new_fig(rect, 4.5, 45, [], [], [], [], "");
% plot_scatter(cs, ms, lw, malpha, cs.blu, log10(xc), log10(nanmean(yc, 1)), ...
%     nanstd(log10(yc), 0, 1));
% plot_scatter(ss, ms, lw, malpha, cs.ora, log10(xs), log10(nanmean(ys, 1)), ...
%     nanstd(log10(ys), 0, 1));
eb = shadedErrorBar(log10(xs), log10(nanmean(ys, 1)), nanstd(log10(ys), 0, 1)); 
eb.patch.FaceColor = cs.ora; eb.mainLine.LineStyle = 'none'; eb.edge(1).LineStyle = 'none'; eb.edge(2).LineStyle = 'none';
eb = shadedErrorBar(log10(xc), log10(nanmean(yc, 1)), nanstd(log10(yc), 0, 1)); 
eb.patch.FaceColor = cs.blu; eb.mainLine.LineStyle = 'none'; eb.edge(1).LineStyle = 'none'; eb.edge(2).LineStyle = 'none';
% plot(log10(xc), log10(nanmean(yc, 1)), 'linewidth', lw, 'color', c.blu); hold on
% plot(log10(xs), log10(nanmean(ys, 1)), 'linewidth', lw, 'color', c.ora);
scatter(log10(xc), log10(nanmean(yc, 1)), ms , 'markerfacecolor', cs.blu,'markeredgealpha',0, 'markerfacealpha', malpha); hold on
scatter(log10(xs), log10(nanmean(ys, 1)), ms , 'markerfacecolor', cs.ora,'markeredgealpha',0, 'markerfacealpha', malpha);
xt = 0:6; xticks(xt);
% xtl = ["10^{3}", "10^4", "10^{5}"]; xticklabels(xtl)
yt = 2:2:6; yticks(yt)
% ytl = ["10^2", "10^3", "10^4", "10^5"]; yticklabels(ytl)
% xlim([3 4.93])
% xlim([3 5])
xlabel("log t (min)")
ylabel("\delta^2 (\mum^2)")
% r = refline(1, -1.05); r.LineStyle = '--'; r.LineWidth = 1.5; r.Color = cs.gra;
% r = refline(2, -3.85); r.LineStyle = '--'; r.LineWidth = 1.5; r.Color = cs.gra;
r = line([1 2], [1.8 3.8]); r.LineWidth = 3; r.LineStyle = '--'; r.Color = 'k';
r = line([2.4 3.4], [2.8 3.8]); r.LineWidth = 3; r.LineStyle = '--'; r.Color = 'k';
ylim([1 6])
xlim([1 4])
setfigprop(800,800,4, 40)

% Msd inset
ms = 150; lw = 3; mlapha = 0.9;
new_fig(rect, 4.5, 45, [], [2.7 5.3], [], [], "");
plot(log10(xc), log10(nanmean(yc, 1)), 'linewidth', lw, 'color', cs.blu); hold on
plot(log10(xs), log10(nanmean(ys, 1)), 'linewidth', lw, 'color', ss.ora);
scatter(log10(xc), log10(nanmean(yc, 1)), ms , 'markerfacecolor', cs.blu,'markeredgealpha',0, 'markerfacealpha', malpha); hold on
scatter(log10(xs), log10(nanmean(ys, 1)), ms , 'markerfacecolor', ss.ora,'markeredgealpha',0, 'markerfacealpha', malpha);
xt = [4 4.5]; xticks(xt)
% xtl = ["10^4", "10^{4.5}"]; xticklabels(xtl)
xtl = ["4", "4.5"]; xticklabels(xtl);
yt = 4; yticks(yt)
ytl = ["10^4"]; yticklabels(ytl)
% hold on; r = refline(alpha1 , log10(myc(t0)) - alpha1 * log10(xc(t0))); r.LineWidth = 2; r.Color = 'k'; r.LineStyle = '-';
% hold on; r = refline(alpha2 , log10(mys(t0)) - alpha2 * log10(xs(t0))); r.LineWidth = 2; r.Color = 'k'; r.LineStyle = '-';
% xlim([log10(xs(t0)) 4.57])
xlim([3.8 4.7])
ylim([3.2 4.3])
% ylim([3.7 4.2])

% Alpha (do this first than Msd inset)
myc = squeeze(nanmean(yc, 1));
mys = squeeze(nanmean(ys, 1));
dt = 1;
numc = log10(myc(1+dt:end)) - log10(myc(1:end-dt));
dnumc = log10(xc(1+dt:end)) - log10(xc(1:end-dt));
alphac = numc ./ dnumc;
nums = log10(mys(1+dt:end)) - log10(mys(1:end-dt));
dnums = log10(xs(1+dt:end)) - log10(xs(1:end-dt));
alphas = nums ./ dnums;

ms = 150;
new_fig([100 100 400 400], 4.5, 45, [], [], [], [], "");
% yyaxis left
new_logxc = (log10(xc(1:end-dt)) + log10(xc(1+dt:end)))/2;
new_logxs = (log10(xs(1:end-dt)) + log10(xs(1+dt:end)))/2;
% plot(new_logxc, alphac, 'x-', 'linewidth', 7, 'color', [cs.colors(1, :) 0.5]);
% plot(new_logxs, alphas, 'x-', 'linewidth', 7, 'color', [cs.colors(2, :) 0.5]);
scatter(new_logxc, alphac, ms , 'markerfacecolor', cs.blu,'markeredgealpha',0, 'markerfacealpha', malpha); hold on
scatter(new_logxs, alphas, ms , 'markerfacecolor', ss.ora,'markeredgealpha',0, 'markerfacealpha', malpha);
% xlim([2.7 4.8])
% ylim([0.8 2])
r = refline(0, 1); r.LineWidth = 3; r.Color = 'k'; r.LineStyle = '--';
ylim([0.7 2])
yticks(0.5:0.5:2)
% yyaxis right
% plot(log10(xc), log10(myc), 'x-', 'linewidth', 4, 'color', cs.colors(1, :));
% plot(log10(xs), log10(mys), 'x-', 'linewidth', 4, 'color', cs.colors(2, :));
% % Linear slope
% t1 = 45; t0 = 26;
% alpha1 = (log10(myc(t1)) - log10(myc(t0))) / (log10(xc(t1)) - log10(xc(t0)));
% hold on; r = refline(alpha1 , log10(myc(t0)) - alpha1 * log10(xc(t0))); r.LineWidth = 1;
% r.Color = 'k';
% t1 = 45; t0 = 26;
% alpha2 = (log10(mys(t1)) - log10(mys(t0))) / (log10(xs(t1)) - log10(xs(t0)));
% hold on; r = refline(alpha2 , log10(mys(t0)) - alpha2 * log10(xs(t0))); r.LineWidth = 1;
% r.Color = 'k';
ylabel('\alpha')
xlabel("log t (min)")
setfigprop(800,800,4, 40)

%% Population image & load position
% Load all cells
load("C:\Users\HyunGyu\mda\CellCluster Class\simulation trajectories\S 2.8 E -65\conf.mat");
% Make periodic
pos_per = pos;
v = diff(pos_per, 1, 1);
for t = 1:size(pos_per,1)-1
    pos_per(t+1,:,:) = pos_per(t,:,:) + v(t,:,:);
    pos_per(t+1,:,:) = pos_per(t+1,:,:) - (pos_per(t+1,:,:) > 300)*300 + (pos_per(t+1,:,:) < 0)*300;
end

% Get neighbor indecies
t0 = 501; N = 990; % Initial time 
pos_a = squeeze(pos_per(t0, :, :));
NN = 6; % Neighbor N
neighbor = nan(990, NN); 
for i = 1:N
    r = pos_a(:, i) - pos_a(:, :); 
    r = squeeze(absvec(r, 1));
    [~, idx] = sort(r);
    neighbor(i, :) = idx(1:NN);
    waitbar(i/N)
end

%% Read conf image and modify
white = [255; 255; 255];
black = [0; 0; 0];
yellow = [255; 255; 0]; % 1
blue = [0; 0; 255]; % 3
green = [0; 100; 0]; % 2
orange = [255; 165; 0]; % 4

img_conf = imread("C:\Users\HyunGyu\mda\CellCluster Class\simulation trajectories\S 2.8 E -65\conf_cpm_S2.8_E-65_seed10101-1.tif", floor(t0/10)+1);
L = length(img_conf);

img_copy = img_conf;
for i = 2:size(img_conf, 1)-1
    for j = 2:size(img_conf, 2)-1
        a = squeeze(img_copy(i, j, :));
        if isequal(a, black)
            n1 = img_copy(i+1,j,:);
            img_conf(i,j,:) = n1;
        end
    end
end

img_copy = img_conf;
for i = 2:size(img_conf, 1)-1
    for j = 2:size(img_conf, 2)-1
        a = squeeze(img_copy(i, j, :));
        if isequal(a, black)
            n1 = img_copy(i,j+1,:);
            img_conf(i,j,:) = n1;
        end
    end
end

img_copy = img_conf;
for i = 2:size(img_conf, 1)-1
    for j = 2:size(img_conf, 2)-1
        a = squeeze(img_copy(i, j, :));
        if isequal(a, black)
            n1 = img_copy(i+1,j+1,:);
            img_conf(i,j,:) = n1;
        end
    end
end

% % Extract all unique colors (from cell centers)
% color_conf = [];
% for i = 1:990
%     x = squeeze(pos_per(t0, 1, i)) / 300 * L; 
%     y = squeeze(pos_per(t0, 2, i)) / 300 * L;
%     x = floor(x) + 1;
%     y = floor(L-y) + 1;
%     col = squeeze(img_conf(x, y, :)); col = col';
%     a = sum(ismember(color_conf, col), 2);
%     if ~any(a==3)
%         color_conf = [color_conf; col];
%     end
% end

% Extract all unique colors (from all positions)
color_conf = [nan nan nan];
for i = 1:L
    for j = 1:L
        col = squeeze(img_conf(i, j, :)); col = col';
%         a = sum(isequal(col, color_conf), 2);
        a = col == color_conf;
        a = sum(a, 2);
        if ~any(a==3)
            color_conf = [color_conf; col];
        end
    end
end

% Count the how many times each color appear
color_ind = zeros(length(color_conf), 1);
for i = 1:L
    for j = 1:L
        col = squeeze(img_conf(i, j, :)); col = col';
        A = col == color_conf;
        A = sum(A, 2);
        idx = find(A == 3);
        color_ind(idx) = color_ind(idx) + 1;
    end
end

[a, idx] = sort(color_ind, 'descend');

color_picked = color_conf(idx(1:10), :);

figure; colorbar; colormap(double(color_picked)/255)

figure; scatter(1:10, 1:10, 300, double(color_picked)/255, 'filled')

% figure; scatter(1:8, 1:8, 300, cs.colors, 'filled')

figure; imshow(img_conf)

%%
% Remove stupid colors
img_copy = img_conf;
for i = 2:size(img_conf, 1)-1
    for j = 2:size(img_conf, 2)-1
        a = squeeze(img_copy(i, j, :)); a = a';
        A = a == color_picked; A = sum(A, 2);
        if ~any(A == 3)
            b = img_conf(i-1:i+1, j-1:j+1, :);
            b = reshape(b, 9, 3);
            b(5,:) = [];
            N = zeros(10,1); % Counter for the colors near target (i, j)
            for n = 1:8
                bn = b(n, :);
                A = sum(bn == color_picked,2); 
                idx = find(A == 3);
                if ~isempty(idx)
                    N(idx) = N(idx) + 1;
                end
            end
            [~, idx] = max(N);
            img_conf(i,j,:) = color_picked(idx, :);
        end
    end
end

img_copy = img_conf;

% Scale the colors (up)
factor = 1;
offset = 30;
img_copy = uint8(min(double(img_copy) * factor + offset, 255));
% img_conf2 = uint8(double(img_conf2) * 0.90);

% Water it down 
a = 0.8;
A = double(img_copy);
mA = nanmean(img_copy, 3);
dA = A - mA;
A = mA + a * dA;
img_copy = uint8(A);

figure; imshow(img_copy)

img_conf = img_copy;



%% Checking
% t0 = 501;
% img = imread("C:\Users\HyunGyu\mda\CellCluster Class\simulation trajectories\S 2.8 E -65\conf_cpm_S2.8_E-65_seed10101-1.tif", floor(t0/10)+1);
img = img_conf;
figure; imshow(img);
hold on 
x = squeeze(pos_per(t0, 1, :)) / 300 * L; % Image is scaled to 500 (actually 300x300)
y = squeeze(pos_per(t0, 2, :)) / 300 * L;
scatter(x, L-y, 10, 'filled');
for i = 1:990
    xx = x(i);
    yy = y(i);
    text(xx, L-yy, num2str(i), 'color', 'r');
end

%% 3D trajectory (confluent) - Tube
target = 399; % Target cell index
Lim = 200 / 3;
lw = 6;
colors = [cs.blu; cs.ora; cs.yel; cs.pur; cs.sky; cs.gra];
index = squeeze(neighbor(target, :));
% Extract the colors from the image itself
col_ext = [];
for i = index
    x = squeeze(pos_per(t0, 1, i));
    y = squeeze(pos_per(t0, 2, i)); y = 300 - y; % 300 is the actual space width (not the image size)
    x = round(x / 300 * L);
    y = round(y / 300 * L);
    col = squeeze(img_conf(y, x, :));
    col_ext = [col col_ext];
end
new_fig([50 50 1000 1000], 4.5, 45, [], [], [], [], ""); hold on 
x = squeeze(pos_per(t0:(t0+450), 1, index));
y = squeeze(pos_per(t0:(t0+450), 2, index)); y = 300 - y;
% z = (0:450) / 30; z = z';
z = (0:450); z = z';
mx = nanmean(x(1,:), 2); 
my = nanmean(y(1,:), 2);

interpN = 4; twfactor = 30; sp = 0.0005;
zfactor = 0.3; twoffset = 0;
% Draw each after smoothing them
for ii = 1:NN
    ind = ~isnan(x(:,ii)) & ~isnan(y(:,ii));
    x_fit = fit(z(ind), x(ind, ii), 'smoothingspline', 'smoothingparam', sp);
    y_fit = fit(z(ind), y(ind, ii), 'smoothingspline', 'smoothingparam', sp);
%     x_smooth = x_fit(z(ind));
%     y_smooth = y_fit(z(ind));
%     plot3(x_smooth, y_smooth, z(ind), 'linewidth', lw, 'color', colors(ii,:))
    z_new = z(ind);
    dz = z_new(2) - z_new(1);
    z_new = z_new(1): dz/interpN : z_new(end); z_new = z_new';
    x_smooth = x_fit(z_new);
    y_smooth = y_fit(z_new);
    z_new = z_new * zfactor;
    vx = diff(x_smooth,1);
    vy = diff(y_smooth,1);
    v = sqrt(vx.^2 + vy.^2);
    tw = v * twfactor + twoffset; tw = [tw; tw(end)];
    tubeplot(x_smooth, y_smooth, z_new, tw, ones(length(z_new),1), 30);
    
end
% Color the tubes
surface = findobj(gcf, 'type', 'surface');
for i = 1:NN
%     ii = NN-i + 1;
%     surface(i).EdgeColor = 'none';
%     surface(i).FaceColor = colors(ii,:);
%     surface(i).FaceAlpha = 0.8;
    surface(i).EdgeColor = 'none';
    surface(i).FaceColor = col_ext(:, i);
    surface(i).FaceAlpha = 0.8;
end
camlight right
camlight left

xlim([mx-Lim mx+Lim])
ylim([my-Lim my+Lim])
% zlim([0 15])

xt = (mx-Lim:Lim/2:mx+Lim); xticks(xt);
yt = (my-Lim:Lim/2:my+Lim); yticks(yt);
% xtl = string(-Lim:Lim/2:Lim); xticklabels(xtl);
% ytl = string(-Lim:Lim/2:Lim); yticklabels(ytl); 
xticklabels(repelem("", 5));
yticklabels(repelem("", 5));
dz = z_new(2) - z_new(1); 
z_new = z_new(1):dz:(z_new(end)+dz*1000);
zt = z_new(1:150*interpN:end);
z_hour = z_new / zfactor / 30;
ztl = string(z_hour(1:150*interpN:end)*60);
zticks(floor(zt));
zticklabels(ztl);
zlim([0 15*30*zfactor])
zlabel('t (min)')

% Overlay bottom
[X, Y] = meshgrid(0:L-1, 0:L-1);
X = X / L * 300; % The original image is stretched to L x L. 300 is actual pixel width (height) of the simulation space
Y = Y / L * 300;
Z = zeros(L, L);
grid on
sur = surf(X, Y, Z, img_conf); sur.EdgeColor = 'none'; sur.FaceColor = 'texturemap';
set(gcf, 'pos', [86,125,871,839])
set(gca, 'view', [-133,18])
setfigprop(1200,1200,4, 40)

%% Single cell image
% Import the image
A = imread("C:\Users\HyunGyu\mda\CellCluster Class\simulation trajectories\S 3.5 E -50\single_cell_img(92,99center).png"); L = length(A);
L = round(2*Lim*2.5); hL = round(L/2);
img_single = 255*ones(L,L,3); img_single = uint8(img_single);
img_single(hL-50:hL+50, hL-50:hL+50, :) = A(99-50:99+50, 92-50:92+50, :);

% Modify image color
factor = 1.5;
for i = 1:size(img_single, 1)
    for j = 1:size(img_single, 2)
        a = squeeze(img_single(i, j, :));
        if ~isequal(a, white) && ~isequal(a, black)
            if isequal(a, yellow) 
                img_single(i, j, :) = min(cs.blu * 255 * factor, 255);
%                 img_single(i, j, :) = 150*ones(3,1);
            else
                img_single(i, j, :) = [50 50 50];
            end
        elseif isequal(a, white)
            img_single(i, j, :) = black;
        end
    end
end

%% 3D trajectory (single)
load("C:\Users\HyunGyu\mda\CellCluster Class\simulation trajectories\S 2.8 E -65\single.mat");
index = 1:5;
t0 = 501;
new_fig([50 50 1000 1000], 4.5, 45, [], [], [], [], ""); hold on 
x = squeeze(pos(t0:(t0+450), 1, index)); x = x - x(1, :);
y = squeeze(pos(t0:(t0+450), 2, index)); y = y - y(1, :);
% z = (0:450) / 30; z = z';
z = 0:450; z = z';

interpN = 4; twfactor = 30; sp = 0.0005;
zfactor = 0.3; twoffset = 0;
% Draw each after smoothing them
for ii = 1:5
    ind = ~isnan(x(:,ii)) & ~isnan(y(:,ii));
    x_fit = fit(z(ind), x(ind, ii), 'smoothingspline', 'smoothingparam', sp);
    y_fit = fit(z(ind), y(ind, ii), 'smoothingspline', 'smoothingparam', sp);
%     x_smooth = x_fit(z(ind));
%     y_smooth = y_fit(z(ind));
%     plot3(x_smooth, y_smooth, z(ind), 'linewidth', lw, 'color', colors(ii,:))
    z_new = z(ind);
    dz = z_new(2) - z_new(1);
    z_new = z_new(1): dz/interpN : z_new(end); z_new = z_new';
    x_smooth = x_fit(z_new);
    y_smooth = y_fit(z_new);
    z_new = z_new * zfactor;
    vx = diff(x_smooth,1);
    vy = diff(y_smooth,1);
    v = sqrt(vx.^2 + vy.^2);
    tw = v * twfactor + twoffset; tw = [tw; tw(end)];
    tubeplot(x_smooth, y_smooth, z_new, tw, ones(length(z_new),1), 30);
    
end
% Color the tubes
surface = findobj(gcf, 'type', 'surface');
for i = 1:5
    ii = 6-i;
    surface(i).EdgeColor = 'none';
    surface(i).FaceColor = colors(ii,:);
    surface(i).FaceAlpha = 0.8;
end
camlight right
camlight left
camlight 

xt = (-Lim:Lim/2:Lim); xticks(xt);
yt = (-Lim:Lim/2:Lim); yticks(yt);
% xtl = string(-Lim:Lim/2:Lim); xticklabels(xtl);
% ytl = string(-Lim:Lim/2:Lim); yticklabels(ytl); 
xticklabels(repelem("", 5));
yticklabels(repelem("", 5));
dz = z_new(2) - z_new(1); 
z_new = z_new(1):dz:(z_new(end)+dz*1000);
zt = z_new(1:150*interpN:end);
z_hour = z_new / zfactor / 30;
ztl = string(z_hour(1:150*interpN:end) * 60);
zticks(floor(zt));
zticklabels(ztl);
zlim([0 15*30*zfactor])
xlim([-Lim Lim])
ylim([-Lim Lim])
zlabel('t (min)')
% Overlay bottom
[X, Y] = meshgrid(-L/2:L/2-1, -L/2:L/2-1);
X = X / 2.5; % The original image is stretched to 250x250
Y = Y / 2.5;
Z = zeros(L, L);
sur = surf(X, Y, Z, img_single); sur.EdgeColor = 'none'; sur.FaceColor = 'texturemap';
grid on
% zlim([0 15])
set(gcf, 'pos', [86,125,871,839])
set(gca, 'view', [131,13])
setfigprop(1200,1200,4, 40)

%%
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

plot3d_traj(cs, nind(1), 1:701, rect, 200, viewangle, cs.colors(1,:));
plot3d_traj(cs, nind(2), 1:701, rect, 200, viewangle, cs.colors(2,:));
plot3d_traj(cs, nind(3), 1:701, rect, 200, viewangle, cs.colors(3,:));
plot3d_traj(cs, nind(4), 1:701, rect, 200, viewangle, cs.colors(4,:));
plot3d_traj(cs, nind(5), 1:701, rect, 200, viewangle, cs.colors(2,:));


%% Save all figures
figHandles = findall(0,'Type','figure'); 

savefig(figHandles, "figures\simul_conf_single_dec_24.fig");






%%
rect = [500 500 700 700];
new_fig(rect, 2.5, 25, [3.5 4.5], [], [], [], "");

% plot(log10(120*(0:c.T-1)), log10(nanmean(msd_c, 1)), 'x'); hold on
% plot(log10(120*(0:s.T-1)), log10(nanmean(msd_s(1:78,:), 1)), 'x');

shadedErrorBar(log10(120*(0:ss.T-1)), log10(nanmean(msd_ss(1:78,:), 1)), nanstd(log10(msd_ss(1:78,:)), 0, 1));
shadedErrorBar(log10(120*(0:cs.T-1)), log10(nanmean(msd_cs, 1)), nanstd(log10(msd_cs), 0, 1));

plot(log10(120*(0:cs.T-1)), log10(nanmean(msd_cs, 1)), 'x'); hold on
plot(log10(120*(0:ss.T-1)), log10(nanmean(msd_ss(1:78,:), 1)), 'x');

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
malpha = 0.8;
colors = parula(20);
myc = squeeze(nanmean(yc, 1));
mys = squeeze(nanmean(ys, 1));
new_fig(rect, 4.5, 45, [], [2.7 5.3], [], [], "");
for dt = 1
numc = log10(myc(1+dt:end)) - log10(myc(1:end-dt));
dnumc = log10(xc(1+dt:end)) - log10(xc(1:end-dt));
alphac = numc ./ dnumc;
nums = log10(mys(1+dt:end)) - log10(mys(1:end-dt));
dnums = log10(xs(1+dt:end)) - log10(xs(1:end-dt));
alphas = nums ./ dnums;
% Confluent
new_logxc = (log10(xc(1:end-dt)) + log10(xc(1+dt:end)))/2;
% plot(new_logxc, alphac, '-', 'linewidth', lw, 'color', cs.blu);
scatter(new_logxc, alphac, ms, 'markerfacecolor', cs.blu, 'markerfacealpha', malpha, 'markeredgealpha', 0);
% Single
new_logxs = (log10(xs(1:end-dt)) + log10(xs(1+dt:end)))/2;
% plot(new_logxs, alphas, '-', 'linewidth', lw, 'color', cs.ora);
scatter(new_logxs, alphas, ms, 'markerfacecolor', cs.ora, 'markerfacealpha', malpha, 'markeredgealpha', 0);

% xlim([2.7 4.8])
% xlim([3.4 5])
xlim([2.9 5.5])
ylim([0.8 2.0])
end
r = refline(0, 1); r.Color = 'k';
ylabel('\alpha')
xlabel('t (sec)')

%% Correlation map (confluent)
% idx1 = 1;
% idx1 = 40;

% meancorr = nan(length(edges), length(edges), 36);
% meancou = nan(length(edges), length(edges), 36);

% for idx = 77:2:121
%     load(['pos_' num2str(idx) '.mat'])



%     load('C:\Users\HyunGyu\mda\CellCluster Class\simulation trajectories\S 2.8 E -65\conf.mat');

% Go to "E:\cpm_morpheus\confluent 990 cells\pos a"
% aa = [7    18    29    40    51    62    73    84    95   106   117];       % E = -65 fixed and sweep S
aa = [1    12    23    34    45    56    67    78    89   100   111];        % E = -5 fixed and sweep S;
aa = aa(1:3:end);
idx = 1;
for aaa = aa

    load(['pos_' num2str(aaa) '.mat'])

%         pos = pos(:, :, 1:end);
    pos = permute(pos, [3,2,1]);
    pos = reshape(pos, [1 size(pos)]);
    cs = CellCluster(pos);

    edges = -41:2:41; 
    celln = randperm(990, 400); % Pick only part of the cells
    Time = 2100:100:8000;
%     Time = 2100:100:2200;
    Mcor = nan(length(Time), length(edges), length(edges));
    Ncou = nan(length(Time), length(edges), length(edges));
    for i = 1:length(Time)
        t = Time(i);
        timespan = t : (t+100);
        [Mccorrel, Ncount] = correl_map_periodic(cs, edges, timespan, celln);
        Ncount(21, 21) = 0;
        Mcor(i, :, :) = Mccorrel;
        Ncou(i, :, :) = Ncount;
%         disp(t)
         waitbar(i/length(Time))
    end

    meancorr(:, :, idx) = squeeze(nanmean(Mcor, 1));
    meancou(:, :, idx) = squeeze(nanmean(Ncou, 1));
    
    disp(idx)
    
    idx = idx + 3;
    
    

end



%%
figure; 
heatmap(Mccorrel', 'gridvisible', 'off');
% 
figure;
heatmap(Ncount', 'gridvisible', 'off');


%% Diamond heatmap (scatter plot)

margin = 16; 
hs = floor((length(edges)+1)/2);

figure('pos', [50 50 500 500]);
% aa = 1:16;
% for i = 1:length(aa)
%     Mccorrel = meancorr(:, :, i);
%     Ncount = meancou(:, :, i);

% Mean correlation
M = Mccorrel;
M(hs,hs) = 0;
% M(Ncount == 0) = 0;
M = M(hs-margin:hs+margin,hs-margin:hs+margin)';
colorlim = 0.4;
M = M / colorlim;
M(M > 1) = 1;
M(M < -1) = -1;

% Count
% M = Ncount;
% M(hs, hs) = 0;
% M = M(hs-margin:hs+margin,hs-margin:hs+margin)';
% maxcount = 12000;
% M = M / maxcount;

% Draw
% figure('pos', [579,466,561,669]); hold on 
edgesA = edges(1:end-1) + (edges(2) - edges(1))/2;
edgesA = edgesA;
edgesA = edgesA(hs-margin:hs+margin);
[x, y] = meshgrid(edgesA, edgesA);
x_A = x(1:end-1, 1:end-1) + (edgesA(2) - edgesA(1))/2;
y_A = y(1:end-1, 1:end-1) + (edgesA(2) - edgesA(1))/2;
% figure('pos', [100 100 401 409]); hold on 
% scatter(x(:), y(:), 100 , M(:), 'filled', 'd')
A = (M(1:end-1, 1:end-1) + M(1:end-1, 2:end) + M(2:end, 1:end-1) + M(2:end, 2:end)) / 4;
scatter(x_A(:), y_A(:), 100 , A(:), 'filled', 'd'); hold on 
scatter(x(:), y(:), 100 , M(:), 'filled', 'd')
colormap(parula)
% colormap(spring)
xt = -80 : 20 : 80;
yt = -80 : 20 : 80;
xticks(xt)
yticks(yt)
set(gca, 'fontsize', 40, 'linewidth', 4)
box on
% xlim([-30 30])
% ylim([-30 30])
xlim([-27 27])
ylim([-27 27])
xlabel('x(\mum)')
ylabel('y(\mum)')
setfigprop(601, 687, 4, 40)


% imwrite(gcf, ['corr_' num2str(in) '.png'])
% saveas(gcf,  ['corr_' num2str(i) '.png'])
% clf

% end


%% Color bar
figure('pos', [100 100 144 539]); 
cb = colorbar;
cb.FontSize = 40;
cb.TickLength = 0.02;
cb.Ticks = 0:0.25:1;
cb.TickLabels = string(-0.4:0.2:0.4);
% cb.TickLabels = string((0:250:1000)/1000);
cb.LineWidth = 3;
cb.Position = [0.0394 0.0558 0.2200 0.9000];
cb.AxisLocation = 'in';
axis off
colormap(parula)
% colormap(spring)

%% Axis (Load the cropped correlation 'map' as an image and draw axis around it)
% img1 = imread("C:\Users\HyunGyu\Desktop\Paper writings\jan 18\correl(edit).png"); % Read saved image (edited in PowerPoint)
% img1 = imread("C:\Users\HyunGyu\Desktop\Paper writings\jan 18\count(edit).png"); % Read saved image (edited in PowerPoint)
% img1 = imread("C:\Users\HyunGyu\Desktop\Paper writings\feb 15\ncount.png");
% img2 = imread("C:\Users\HyunGyu\Desktop\Paper writings\feb 15\mcorrel.png");

% img_single = imread("C:\Users\Joseph\Desktop\Paper Figures\dec 27\img2(exp).png");

img1 = imread("E:\cpm_morpheus(heterogeneous)\Chain formation analysis\May 10 Figures\mcorrel.png");
img2 = imread("E:\cpm_morpheus(heterogeneous)\Chain formation analysis\May 10 Figures\mcount.png");

H = length(img1) - 1;

fs = 90; lw = 6;
figure;
% figure('pos', [100 100 500 500]);
imshow(img1);
axis on
% set(gca, 'fontsize', 20, 'linewidth', 2)
xt = 1:(H-1)/4:H; xticks(xt);
xtl = string(-80:40:80); xticklabels(xtl);
yt = xt; yticks(yt);
ytl = string(80:-40:-80); yticklabels(ytl);
% xlabel('x(\mum)')
% ylabel('y(\mum)')
set(gca, 'fontsize', fs, 'linewidth', lw)
set(gcf, 'pos', [197,42,1564,1313])
% set(gca, 'outerposition', [0.01,-0.01,0.93,1.047])
% setfigprop(800,657, 4, 40)
% set(gca, 'fontsize', fs, 'linewidth', lw)
% set(gcf, 'pos', [197,42,1564,1313])
% % set(gca, 'outerposition', [0.08,0.13,0.85,0.89]) % To avoid labels being cut when in full screen
% % set(gcf, 'pos', [197,42,1564,1313])


H = length(img2);
figure;
% imshow(img_single);
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
% set(gca, 'outerposition', [0.01,-0.01,0.93,1.047])
% set(gca, 'outerposition', [0.08,0.13,0.85,0.89]) % To avoid labels being cut when in full screen



%% Persistence time for each (200)
Coeff = [];
% figure;
for i = 1:200
    y = auto_correl_cs(i, :);
    y = y(1:100);
    x = 0:100-1; 
    x = x * 2;
    f = fit(x', y', 'exp2', 'lower', [0 -3 0 -3], 'upper', [1 -0.01 1 -0.01]);
%     plot(f, x, y)
%     pause(0.1)
%     clf
    coeff = coeffvalues(f);
    Coeff = [Coeff; coeff];
    waitbar(i / 200);
end
a = - 1 ./ Coeff(:, [2 4]);
a1 = max(a, [], 2);
tau_c = a1; % Minutes

Coeff = [];
% figure;
for i = 1:200
    y = auto_correl_ss(i, :);
    y = y(1:100);
    x = 0:100-1; 
    x = x * 2;
    f = fit(x', y', 'exp2', 'lower', [0 -3 0 -3], 'upper', [1 -0.01 1 -0.01]);
%     plot(f, x, y)
%     pause(0.1)
%     clf
    coeff = coeffvalues(f);
    Coeff = [Coeff; coeff];
    waitbar(i / 200);
end
a = - 1 ./ Coeff(:, [2 4]);
a1 = max(a, [], 2);
tau_s = a1; % Minutes

% Exclude outliers
% tau(tau > 100) = [];

%%
% Exclude 
Coeff_copy = Coeff;
Coeff = [];
for i = 1:200
    if Coeff_copy(i, 1) > 0.2
        Coeff = [Coeff; Coeff_copy(i, :)];
    end
end

% Sort the coeffcients
for i = 1:200
    A = Coeff(i, 1);
    B = Coeff(i, 3);
    t1 = Coeff(i, 2);
    t2 = Coeff(i, 4);
    if A < B
        Coeff(i, 1:2) = Coeff(i, 3:4);
        Coeff(i, 3:4) = [A t1];
    end
end

% ptime_cs = - 1./ Coeff(:, 2);
% ptime_ss = - 1 ./ Coeff(:, 2);

ptime_cs = - 1./ Coeff(:, 2);
ptime_ss = - 1./ Coeff(:, 2);

%% Speed, persistent length
v_c = abs(diff(cs.pos, 1, 4)) * cs.scale /2; % um / min
v_s = abs(diff(ss.pos, 1, 4)) * ss.scale /2; % um / min

v_c = absvec(v_c, 3);
v_s = absvec(v_s, 3);

mv_c = squeeze(nanmean(v_c, 4))';
mv_s = squeeze(nanmean(v_s, 4))';

plength_c = tau_c .* mv_c;
plength_s = tau_s .* mv_s;


%% Diffusion coefficients
t = 0:5000-1;
t = t * 2; % In minites
logt = log10(t);
% figure; plot(logt, log10(msd_cs)); 
figure; plot(logt, log10(msd_ss))
ylim([0 7]); xlim([0 4])

% y = log10(msd_cs);
y = log10(msd_ss);

dy = y(:, 3) - y(:, 2);
dx = logt(3) - logt(2);
a = dy ./ dx;
y0 = y(:, 2);
x0 = logt(2);
logD = y0 - a .* x0;
D = 10.^logD;













