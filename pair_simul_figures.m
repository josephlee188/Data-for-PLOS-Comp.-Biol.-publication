% s = 8; e = 11;
trange = 1000:4000; % Using whole range could be too heavy
pos_pair_simul = permute(pos(trange, :, :, :), [3 4 2 1]);
% pos_pair_simul = pos_pair_simul / 10 * 30; % Scale to match um
ps = CellCluster(pos_pair_simul);
% ps.scale = 1;
ps.scale = 1;

dth = ang_vel(ps);
D = cell_dist(ps);
th = ang(ps);

dt = 5;
dth_2 = nan(2, ps.N, ps.T-dt); % dth_2 with new interval
for t = 1:ps.T-dt
    dth_2(:, :, t) = th(:, :, t+dt) - th(:, :, t);
end
dth_2 = dth_2 / dt;

%% 3D trajectory 
cell_ind = 2; 
% dur = 450; t0 = 2100; t0 = 1380;
dur = 3000; t0 = 1;
timespan = t0 : (t0+dur);
rect = [1400 200 700 850];
viewangle = [327.6,379.0,41.9];
figure('pos', rect); hold on 
% lim = 50;
lim =30;
plot3d_traj(ps, cell_ind, timespan, rect, lim, viewangle, []);
% xt = -100:50:100; xticks(xt); xtl = ["-100" "" "0" "" "100"]; xticklabels(xtl);
% yt = -100:50:100; yticks(yt); ytl = ["-100" "" "0" "" "100"]; yticklabels(ytl);
zlabel('t')
xt = -30:15:30; xticks(xt);
yt = -30:15:30; yticks(yt);
% img = imread('..\COMBINED_pair_Oct8_Oct23_Oct15(59).tif', 1);
% hold on; overlay_image(p, img, cell_ind, lim, 0, 4);
img = squeeze(IMG(1, :, :, :));
% img = permute(img, [2 1 3]);
img = img(end:-1:1, :, :);
[X,Y] = meshgrid(-49:50, -49:50); Z = zeros(100, 100); 
X = X / 2; Y = Y / 2;
X = X + 2;
Y = Y + 0;
s = surf(X, Y, Z, img); s.EdgeColor = 'none';
xlabel('x')
ylabel('y')

%% 3D Tube
scale = 2.6;

cell_ind = 2; 
dur = 450; t0 = 1380;
timespan = t0 : (t0+dur);
rect = [0 0 700 850];
lim = 13;
viewangle = [183.4,124.5,72.9];
new_fig(rect, 4.5, 45, [], [], [], [], "");
interpN = 4; 
% twfactor = 10; 
twfactor = 7;
zfactor = 0.08;
sp = 0.008;
tube_surf(ps, cell_ind, timespan, rect, lim, viewangle, twfactor, interpN, sp, zfactor);
% Overlay bottom
img_ground = zeros(200,200,3); Z0 = zeros(200,200); 
[X,Y] = meshgrid(-99:100, -99:100);
s = surf(X,Y,Z0, img_ground); s.EdgeColor = 'none';
img = squeeze(IMG(1, :, :, :));
img = img(end:-1:1, :, :);
[X,Y] = meshgrid(-49:50, -49:50); Z = zeros(100, 100)+0.01; 
X = X / scale; Y = Y / scale;
X = X + 0;
Y = Y + 0;
s = surf(X, Y, Z, img); s.EdgeColor = 'none';
% xlabel('x')
% ylabel('y')
xt = -20/3:20/3:20/3; xticks(xt)
yt = -20/3:20/3:20/3; yticks(yt)
xticklabels(repmat("", [3 1]))
yticklabels(repmat("", [3 1]))
zlim([0 15*30*zfactor])

camroll(-90)
setfigprop(1700,1200,4,40)
set(gca, 'innerpos', [0.13,0.11,0.78,0.73])
ztickangle(-90)

%% Pairs
dth_2_smooth = smoothing(ps, dth_2, 0.85);
D_smooth = smoothing(ps, D, 0.1);
th_smooth = smoothing(ps, th, 0.75);

% theta vs. time
time = (0:ps.T-1)/30;
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], ""); hold on 
% yyaxis left; 
ylabel('\theta (rad)')
tspan = timespan - timespan(1); tspan = tspan * 2;
plot(tspan, squeeze(th_smooth(1, cell_ind, timespan) - th_smooth(1, cell_ind, timespan(1)) + pi/6), 'linewidth', 6, 'color', ps.colors(1,:));
% plot(tspan, squeeze(th(1, cell_ind, timespan)), 'linewidth', 6, 'color', ps.colors(1,:));
yt = -4*pi : 2*pi : 4*pi; yticks(yt);
ytl = string(-4:2:4) + "\pi"; ytl(3) = '0'; yticklabels(ytl);
% r = refline(0,0); r.LineWidth = 2; r.LineStyle = ':';
ylim([-5*pi 3*pi+1])
xlim([0 900])
xt = 0:300:900; xticks(xt)
% ytl = {'-6\pi', '0', '6\pi'}; yticklabels(ytl);
% xlabel('t (min)')
box on
scatter(198, 6, 500, '^', 'k', 'filled')
scatter(784, -11, 500, 'v', 'k', 'filled')
r = refline(0,0); r.LineWidth = 3; r.LineStyle = '--'; r.Color = p.gra;
setfigprop(1000, 450, 4, 40)

% omega vs. D
ms = 100; lw = 4;
% ms = 60; lw = 3;
time = (0:ps.T-1)/30;
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], "");
yyaxis left; ylabel('|\omega| (rad/min)')
plot(tspan, squeeze(abs(dth_2_smooth(1, cell_ind, timespan)) * 30 / 60), 'linewidth', lw,  'color', [ps.gra 0.9]);
scatter(tspan, squeeze(abs(dth_2(1, cell_ind, timespan)) * 30 / 60), ms,  'markerfacecolor', ps.gra, ...
    'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', ps.gra)
yt = 0 : 0.04 : 0.08; yticks(yt);
% ytl = {'0', '\pi', '2\pi', '3\pi', '4\pi'}; yticklabels(ytl);
ylim([0 0.09])
yyaxis right; ylabel('\itd (\mum)')
plot(tspan, squeeze(D_smooth(1, cell_ind, timespan)) * 3, 'linewidth', lw, 'color', ps.gre);
scatter(tspan, squeeze(D(1, cell_ind, timespan)) * 3, ms , 'markerfacecolor', ps.gre, ...
        'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', ps.gre)
% xlabel('t (min)')
ylim([8 14])
yt = 8:2:14; yticks(yt);
xt = 0:300:900; xticks(xt)
xlim([0 900])
box on
setfigprop(1040, 450, 4, 40)
set(gca, 'innerpos', [0.15,0.2502,0.711858981952835,0.749814818957713])

% omega vs. D (inset)
ms = 200; lw = 8;
% ms = 60; lw = 3;
time = (0:ps.T-1)/30;
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], "");
yyaxis left; 
% ylabel('|\omega| (rad/min)')
plot(tspan, squeeze(abs(dth_2_smooth(1, cell_ind, timespan)) * 30 / 60), 'linewidth', lw,  'color', [ps.gra 0.9]);
scatter(tspan, squeeze(abs(dth_2(1, cell_ind, timespan)) * 30 / 60), ms,  'markerfacecolor', ps.gra, ...
    'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', ps.gra)
yt = 0 : 0.04 : 0.08; yticks(yt);
% ytl = {'0', '\pi', '2\pi', '3\pi', '4\pi'}; yticklabels(ytl);
ylim([0 0.09])
yyaxis right; 
ylabel('\itd (\mum)')
plot(tspan, squeeze(D_smooth(1, cell_ind, timespan)) * 3, 'linewidth', lw, 'color', ps.gre);
scatter(tspan, squeeze(D(1, cell_ind, timespan)) * 3, ms , 'markerfacecolor', ps.gre, ...
        'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', ps.gre)
% xlabel('t (min)')
ylim([8 14])
yt = 8:2:12; yticks(yt);
xt = 300:150:600; xticks(xt)
xlim([300 600])
box off
setfigprop(1040, 400, 4, 40)


%%
rect = [50 50 600 500];
new_fig(rect, 2.5, 35, [0 6], [], 0:2:6, [], "D period"); histogram(periodD(1,:));
new_fig(rect, 2.5, 35, [0 6], [], 0:2:6, [], "\omega period"); histogram(periodDth(1,:));
new_fig(rect, 2.5, 35, [], [], [], [], "\rho"); histogram(pcorrel(1:57), 'binwidth', 0.075);

% Draw 3D trajectory (with overlay)
cell_ind = 31; 
timespan = 1:668;
rect = [1400 200 700 850];
lim = 30;
viewangle = [327.6,379.0,41.9];
figure('pos', rect); hold on 
plot3d_traj(ps, cell_ind, timespan, rect, lim, viewangle, []);
img = imread('..\COMBINED_pair_Oct8_Oct23_Oct15(59).tif', 1);
hold on; overlay_image(ps, img, cell_ind, lim, 0, 4);

%% Statistics
pc = pair_correlation(ps);
mpc = nanmean(pc, 2);

newD = subtract_smoothing(ps, D, 3e-6);
newDth = subtract_smoothing(ps, dth, 3e-6) * 30; % rad per hour

periodD = get_period_fft(ps, newD); % In hour
periodDth = get_period_fft(ps, newDth); % hour

pcorrel = pearson_correlation(ps, newD, abs(newDth));

% Histograms
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(periodD*60, 'binwidth', 6, 'normalization', 'probability');
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0.0; h.LineWidth = 1.5;
xlabel('\tau_{\itd} (min)')
xticks(0:120:360)
xlim([0 360])
scatter(nanmean(periodD(:)*60), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 470, 4, 40)

new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(pcorrel, 'binwidth', 0.01, 'normalization', 'probability');
h.FaceColor = [0.00,0.45,0.74]; h.EdgeAlpha = 0; h.LineWidth = 1.5;
xlabel('\rho')
xt = [-0.3 0 0.3]; xticks(xt);
xlim([-0.35 0.35])
scatter(nanmean(pcorrel(:)), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 440, 4, 40)

% Histogram of pair correlation
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(mpc, 'binwidth', 0.01, 'normalization', 'probability');
h.FaceColor = [0.85,0.33,0.10]; h.EdgeAlpha = 0; h.LineWidth = 1.5;
xlabel('c_{pair}')
xt = [-0.4 0 0.4]; xticks(xt);
yticks([0, 0.1])
xlim([-0.5 0.5])
scatter(nanmean(mpc(:)), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 470, 4, 40)


%%
figHandles = findall(0,'Type','figure'); 

savefig(figHandles, "figs_pair_dec23.fig");


%% Scatter (save)
figure;
L = 20; i = 1;
for t = t0 : (t0+dur)
% for t = t0
    x = ps.pos(:, cell_ind, 1, t);
    y = ps.pos(:, cell_ind, 2, t);
    mx = nanmean(x); my = nanmean(y);
    scatter(x, y, 50, 'markerfacecolor', 'k')
    for j = 1:2
        text(x(j), y(j), num2str(j), 'fontsize', 20, 'color', 'r')
    end
    xlim([mx-L mx+L])
    ylim([my-L my+L])
    saveas(gca, ['img_' num2str(i) '.png'])
    i = i + 1;
end

%% Load and adjust CPM image sequence
% Read the image stack (original)
IMG = [];
for i = 1:45
    img = imread("C:\Users\HyunGyu\mda\CellCluster Class\simulation trajectories\S 2.8 E -65\pair_cpm_S2.8_E-65_seed10102.tif", i);
    img = reshape(img, [1 size(img)]);
    IMG = cat(1, IMG, img);
end

% Read the modified images
IMG = [];
for i = 1:45
    img = imread(['img_' num2str(i) '.png']);
    img = reshape(img, [1 size(img)]);
    IMG = cat(1, IMG, img);
end

%% Remove boundary (with function)
% Fresh import IMG from .tif
wb = waitbar(0);
for t = 1:45
    img = squeeze(IMG(t, :, :, :));
    IMG(t, :, :, :) = remove_boundary(ps, img, 1);
    waitbar(t/45)
end
close(wb)

%%
% figure; imshow(squeeze(IMG(109, :, :, :)));
white = [255; 255; 255];
black = [0; 0; 0];
yellow = [255; 255; 0]; % 1
blue = [0; 0; 255]; % 3
green = [0; 100; 0]; % 2
orange = [255; 165; 0]; % 4

% Change color to match (blu, ora, yel, pur)
factor = 1.5;
for t = 1:size(IMG,1)
for i = 1:size(IMG, 2)
    for j = 1:size(IMG, 3)
        a = squeeze(IMG(t, i, j, :));
        if ~isequal(a, white) && ~isequal(a, black)
            if isequal(a, yellow) 
                IMG(t, i, j, :) = min(ps.blu * 255 * factor, 255);
            elseif isequal(a, green) 
                IMG(t, i, j, :) = min(ps.ora * 255 * factor, 255);
            else
                IMG(t, i, j, :) = [50 50 50];
            end
        elseif isequal(a, white)
            IMG(t, i, j, :) = black;
        end
    end
end
end

% Save the image
for t = 1:45
%     imshow(squeeze(IMG(t, :, :, :)));
%     saveas(gca,['img_' num2str(t) '.png'])
    imwrite(squeeze(IMG(t, :, :, :)), ['img_' num2str(t) '.png']); % Save matrix as an image directly
end
    
%% One snapshot (beside 3D tube)
img = squeeze(IMG(1, :, :, :));
factor = 1.5;

blu = uint8(min(ps.blu * 255 * factor, 255));
ora = uint8(min(ps.ora * 255 * factor, 255));

for i = 2:size(img, 1)-1
    for j = 2:size(img, 2)-1
        a = squeeze(img(i, j, :));
        if ~isequal(a, white) && ~isequal(a, black)
            if isequal(a, yellow)
                img(i, j, :) = blu;
            elseif isequal(a, green)
                img(i, j, :) = ora;
            end
        end
    end
end

img_copy = img;
for i = 2:size(img, 1)-1
    for j = 2:size(img, 2)-1
        a = squeeze(img_copy(i, j, :));
        if isequal(a, black)
            n1 = img_copy(i+1,j,:);
            %             n2 = img(i-1,j,:);
            %             n3 = img(i,j+1,:);
            %             n4 = img(i,j-1,:);
            img(i,j,:) = n1;
        end
    end
end

img_copy = img;
for i = 2:size(img, 1)-1
    for j = 2:size(img, 2)-1
        a = squeeze(img_copy(i, j, :));
        if isequal(a, black)
            n1 = img_copy(i,j+1,:);
            %             n2 = img(i-1,j,:);
            %             n3 = img(i,j+1,:);
            %             n4 = img(i,j-1,:);
            img(i,j,:) = n1;
        end
    end
end

% Remove stupid colors 
img_copy = img;
for i = 2:size(img, 1)-1
    for j = 2:size(img, 2)-1
        a = squeeze(img_copy(i, j, :)); a = a';
        if ~isequal(a, white') && ~isequal(a, blu) && ~isequal(a, ora)
            b = img(i-1:i+1, j-1:j+1, :);
            b = reshape(b, 9, 3);
            b(5,:) = [];
            
            wt = sum(b == white', 2);
            bl = sum(b == blu, 2);
            or = sum(b == ora, 2);
            
            wt_n = sum(wt == 3);
            bl_n = sum(bl == 3);
            or_n = sum(or == 3);
            
            [~, ind] = max([wt_n, bl_n, or_n]);
            
            if ind == 1
                img(i,j,:) = white;
            elseif ind == 2
                img(i,j,:) = blu;
            elseif ind == 3
                img(i,j,:) = ora;
            end
%             b = img(i+1,j,:);
%             c = img(i-1,j,:);
%             d = img(i,j+1,:);
%             e = img(i,j-1,:);
%             f = img(i-1,j-1,:);
%             g = img(i+1,j+1,:);
%             h = 

        end
    end
end

figure; imshow(img)



%% Overlay trajectory
ms = 120; lw = 8; malpha = 0.8;
scale = 2.0;
shiftx = 52; shifty = 50;
for dt = [0 30 60 90]
figure;
t00 = t0 + 120; % 4 hour
t = t00 + dt;
T = t - t0 + 1;
imshow(squeeze(IMG(T, :, :, :)));
hold on 
mx = ps.pos(:, cell_ind, 1, t00:t); mx = squeeze(nanmean(mx, 1));
my = ps.pos(:, cell_ind, 2, t00:t); my = squeeze(nanmean(my, 1));
for i = 1:4
x = ps.pos(i, cell_ind, 1, t00:t); x = squeeze(x); x = x - mx(1);
y = ps.pos(i, cell_ind, 2, t00:t); y = squeeze(y); y = y - my(1);
x = x * scale;
y = y * scale; y = - y;
x = x + shiftx;
y = y + shifty;
plot(x, y, 'linewidth', lw, 'color', [ps.colors(i, :) malpha]);
scatter(x(end), y(end), ms, 'markerfacecolor', ps.colors(i, :), 'markeredgecolor', 'k');
end
end

%% Get the mean position from colors
x = nan(4, 100); y = nan(4, 100);
n = ones(4, 1);
for i = 1:100
    for j = 1:100
        a = squeeze(IMG(T, i, j, :)); a = a';
        if isequal(a, uint8(ps.blu * 255))
            x(1, n(1)) = j;
            y(1, n(1)) = i;
            n(1) = n(1) + 1;
        elseif isequal(a, uint8(ps.ora * 255))
            x(2, n(2)) = j;
            y(2, n(2)) = i;
            n(2) = n(2) + 1;
        elseif isequal(a, uint8(ps.yel * 255))
            x(3, n(3)) = j;
            y(3, n(3)) = i;
            n(3) = n(3) + 1;
        elseif isequal(a, uint8(ps.pur * 255))
            x(4, n(4)) = j;
            y(4, n(4)) = i;
            n(4) = n(4) + 1;
        end
    end
end
mx = nanmean(x, 2);
my = nanmean(y, 2);

figure; 
imshow(squeeze(IMG(T, :, :, :)));
hold on 
scatter(mx(4), my(4), 100)






































%% Flatten the theta ???
a = squeeze(th(1, :, :));
% figure('pos', [293,609,693,479]); hold on 
figure('pos', [100,100,694,600]); hold on 
lw = 1.5; alpha = 0.5;
for i = 1:100
    x = (0:450-1) / 30;
    y = a(i, 1:450);
    col = rand(1,3);
    plot(x, y, 'linewidth', lw, 'color', [col alpha]); 
end
set(gca, 'linewidth', 3.7, 'fontsize', 37)
yt = (-12:6:12)*pi; yticks(yt);
ytl = string(-12:6:12) + "\pi"; ytl(3) = "0"; yticklabels(ytl);
ylim([-12*pi 12*pi])
xlim([0 15])
xlabel('time')
ylabel('\theta (rad)')

% Get 'A'
Dt = 1:10:401;
A = nan(length(Dt), ps.N, ps.T);
j = 1;
for dt = Dt
    a = squeeze(th(1, :, 1:dt:end));
    b = diff(a, 1, 2);
    b = abs(b / dt);
%     b = b / dt;
    A(j, :, 1:size(b,2)) = b;
    j = j + 1;
end

a = A(3, :); a(isnan(a)) = [];
a = 2*pi ./ a; % Period
a(a > 450) = nan;
a = a / 30; % In 'Hour' 
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(a, 'binwidth', 0.5, 'normalization', 'probability');
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0;
xlabel('\it{T}')
% xlim([0 0.2])
xlim([0 15])

% Draw histograms of velocity of different intervals
figure; hold on 
col = parula(length(Dt));
for i = 1:15
    a = A(i,:);
    a(isnan(a)) = [];
%     histogram(a, 'EdgeColor', 'none', 'normalization', 'probability', 'binwidth', 0.002);
    [h, edges] = histcounts(a, 'binwidth', 0.004,'normalization', 'probability');
    edges = edges(1:end-1) + edges(1)/2;
    plot(edges, h, 'linewidth', 2, 'color', [col(i,:) 0.2]);
end

figure; hold on 
col = parula(10);
for i = 1:10
    a = A(i,:);
    a(isnan(a)) = [];
    histogram(a, 'EdgeColor', 'none', 'normalization', 'probability', 'binwidth', 0.003, 'facealpha', 0.3, 'facecolor', col(i,:));
end
% xlim([-0.25 0.25])
colorbar

% Make msd th
msd_th = nan(ps.N, ps.T);
for dt = 0:ps.T-1
    t0 = 1:ps.T-dt;
    t1 = t0 + dt;
    x = (th(1, :, t1) - th(1, :, t0)).^2;
    x = squeeze(nanmean(x, 3));
    msd_th(:, dt+1) = x;
    waitbar(dt/ps.T)
end

% Draw msd th
figure; plot(log10(0:ps.T-1), log10(msd_th)); hold on 
plot(log10(0:ps.T-1), log10(nanmean(msd_th,1)), 'linewidth', 3);
refline(2, 0); refline(1, 0);





%% Find the locations of reverse turn 
th_smooth = nan(ps.N, ps.T-1);
for i = 1:ps.N
    y = squeeze(th(1, i, :)); y = y(1:end-1);
    x = 1:length(y);
    f = fit(x', y, 'smoothingspline', 'smoothingparam', 2e-5);
    th_smooth(i, :) = f(x);
    waitbar(i/ps.N)
end

figure; 
plot(0:3000-1, permute(th_smooth, [2 1]));
xlim([0 450])

figure; 
for i = 1:100
    subplot(10, 10, i)
    x = 1:size(th_smooth,2);
    y = th_smooth(i,:);
    plot(x, y); hold on 
    scatter(x, squeeze(th(1, i, 1:end-1)), 1);
end

th_xx = nan(ps.N, ps.T-3);
for i = 1:100
    y = squeeze(th_smooth(i, :));
    ddy = diff(y, 2, 2);
    th_xx(i, :) = ddy;
end

th_x = nan(ps.N, ps.T-2);
for i = 1:100
    y = squeeze(th_smooth(i, :));
    ddy = diff(y, 1, 2);
    th_x(i, :) = ddy;
end

figure; 
for i = 1:100
    subplot(10, 10, i)
%     x = 1:size(th_xx,2);
%     y = th_xx(i,:);
    x = 1:size(th_x,2);
    y = abs(th_x(i,:));
    plot(x, y); 
%     ylim([-5e-3 5e-3])
    yyaxis right
    plot(th_smooth(i, :));
    hold on 
    scatter(1:ps.T, squeeze(th(1, i, :)), 1);
    xlim([0 450])
end

%% Unwrapping the theta vs. time
th = squeeze(th(1, :, :));
T = size(th, 2);
t = 0:T-1; t = t * 2;

%%
t = 0:T-1;
dt = 100;
dth = abs(th(:, dt+1:T) - th(:, 1:T-dt))/dt;

threshold = 0.03;

steep = (dth >= threshold);
th_lower = th;
th_lower(steep) = nan;

figure; 
i = 2;
subplot(2,1,1)
plot(t, th(i, :), 'linewidth', 2)
hold on
plot(t, th_lower(i,:), 'x', 'linewidth', 2)
subplot(2,1,2)
plot(t(1:end-dt), dth(i, :), 'linewidth', 2)
hold on 
l = line([0 t(end)], [threshold threshold]); l.LineWidth = 2; l.LineStyle = '--';

        
 %% Find times of reverse-turn
 Thresh = 0.2;
 t_turn = nan(100,200); 
 t_turn(1, :) = 1;
 for i = 1:200
     y = th(i, :);
     
     N = 1;
     t0 = t_turn(1, i);
%      wb = waitbar(0);
     while t0 < T
         
         n = 1;
         dy = nan(T - t0 + 1, 1);
         for t = t0:T
             y_end = th(i, t);
             y_init = th(i, t0);
             y_straight = (y_end - y_init) / (t - t0 + 0.0001) * ((t0:t) - t0) + y_init;
             dy(n) = nanmean((th(i, t0:t) - y_straight).^2);
             n = n + 1;
         end
         ind = find(dy < Thresh);
         t_new =  ind(end) + t0;
         N = N + 1;
         t_turn(N, i) = t_new;
         t0 = t_new;
%          waitbar(t0/T)
     end
%      close(wb)
    
     waitbar(i/200)
 end
 
%  t_turn(isnan(t_turn)) = [];
%%
figure; 
i = 200;
time = 1:T;
tt = t_turn(:, i); tt(isnan(tt)) = [];
plot(time, th(i, :), 'linewidth', 2); 
hold on
plot(time(tt), th(i, tt), 'x', 'linewidth', 2)
%% Get the slopes
w = nan(100, 200);
for i = 1:200
    a = t_turn(:, i);
    a = a(~isnan(a));
    n = 1;
    for j = 1:length(a)-1
        t0 = a(j);
        t1 = a(j+1);
        y0 = th(i, t0);
        y1 = th(i, t1);
        sl = (y1 - y0) / (t1 - t0);
        w(n, i) = sl;
        n = n + 1;
    end
end

% period = 2*pi ./ w;
%% Histograms
mw = nanmean(abs(w), 1);
figure; histogram(mw)
xlabel('<|\omega|>')

mperiod = 2*pi ./ mw;
figure; histogram(mperiod*2)
xlabel('2\pi/|\omega|')        
xlim([0 720])

new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(mperiod*2, 'binwidth', 6, 'normalization', 'probability');
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0.0; h.LineWidth = 1.5;
xlabel('{\itT} (min)')
xticks(0:240:720)
xlim([0 720])
scatter(nanmean(mperiod*2), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 470, 4, 40)

%% Unwrapping the whole theta (based on initial slope)
th_unwrap = th;
for i = 1:200
    sig = sign(w(1, i));
    a = t_turn(:, i);
    a = a(~isnan(a));
    for j = 2:length(a)-1
        t0 = a(j);
        t1 = a(j+1);
        y = th(i, t0:t1);
        if sig ~= sign(w(j, i))
            th_unwrap(i, t0:t1) = - (y - y(1)) + th_unwrap(i, t0);
        else
            th_unwrap(i, t0:t1) = y - y(1) + th_unwrap(i, t0);
        end
    end
end

%% Plot all unwrapped data
figure('pos', [100,100,694,600]); hold on 
for i=1:200
    plot((0:T-1)*2, th_unwrap(i, :), 'linewidth', 2.4, 'color', [rand(1,3) 0.6]);
end
r = refline(0, 0); r.Color = 'k'; r.LineWidth = 3; r.LineStyle = '--';
xlim([0 1200])
xlabel('t (min)')
ylabel('\theta (rad)')
ylim([-16*pi 16*pi])
yt = [-16*pi : 8*pi : 16*pi]; yticks(yt);
ytl = {'-16\pi' '-8\pi' '0' '8\pi' '16\pi'}; yticklabels(ytl)
xt = 0:300:1200; xticks(xt)
box off
setfigprop(1200, 880, 4, 40)
%%










        
     










