%% Quads
pos_quad_temp = pos_quad;
pos_quad_temp(1, :, :, :) = pos_quad(4, :, :, :);
pos_quad_temp(4, :, :, :) = pos_quad(1, :, :, :);
q = CellCluster(pos_quad_temp(:, :, :, 150:558));

dth = ang_vel(q);
D = cell_dist(q);
th = ang(q);
qc = quad_correlation(q);
dt = 5;
dth_2 = nan(5, 40, q.T-dt);
for t = 1:q.T-dt
    dth_2(:, :, t) = th(:, :, t+dt) - th(:, :, t);
end
dth_2 = dth_2 / dt;

newD = subtract_smoothing(q, D, 3e-6);
newDth = subtract_smoothing(q, dth, 3e-6) * 30; % rad per hour

periodD = get_period_fft(q, newD); % In hour
periodDth = get_period_fft(q, newDth); % hour

pcorrel = pearson_correlation(q, newD, abs(newDth));

rect = [50 50 600 500];
new_fig(rect, 2.5, 35, [0 6], [], 0:2:6, [], "D period"); histogram(periodD);
new_fig(rect, 2.5, 35, [0 6], [], 0:2:6, [], "\omega period"); histogram(periodDth);
new_fig(rect, 2.5, 35, [], [], [], [], "\rho"); histogram(pcorrel, 'binwidth', 0.075);

% Draw 3D trajectory (with overlay)
% cell_ind = 30; 
cell_ind = 11;
timespan = 1:409;
rect = [1400 200 700 850];
lim = 50;
viewangle = [327.6,379.0,41.9];
% Adjustments
q_copy = q;
q_copy.pos(3, cell_ind, 2, :) = q_copy.pos(3, cell_ind, 2, :) - 2.5;
q_copy.pos(2, cell_ind, 2, :) = q_copy.pos(2, cell_ind, 2, :) - 1;
q_copy.pos(2, cell_ind, 1, :) = q_copy.pos(2, cell_ind, 1, :) + 1.5;
q_copy.pos(4, cell_ind, 2, :) = q_copy.pos(4, cell_ind, 2, :) - 0.5;
q_copy.pos(1, cell_ind, 1, :) = q_copy.pos(1, cell_ind, 1, :) + 1.5;
q_copy.pos(1, cell_ind, 2, :) = q_copy.pos(1, cell_ind, 2, :) - 1.5;
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], ""); hold on 
plot3d_traj(q_copy, cell_ind, timespan, rect, lim, viewangle, []);
 img = imread('..\combined_stacks_quad_Oct8_Oct23_Oct15(38)(2).tif', 150);
hold on; overlay_image(q_copy, img, cell_ind, lim, 3, -5);
xt = [-40 -20 0 20 40]; xticks(xt)
% xtl = ["-40" "" "0" "" "40"]; xticklabels(xtl)
xtl = ["" "" "" "" ""]; xticklabels(xtl)
yt = xt; yticks(yt)
ytl = xtl; yticklabels(ytl)
zt = 0:2.5:10.0; zticks(zt);
ztl = ["0" "" "5" "" "10"]; zticklabels(ztl);
% f = fill3([-lim -lim -lim -lim], [-50 50 50 -50], [5 5 15 15], [0.93 0.69 0.13 0.3]);
% f.FaceColor = 'y'; f.FaceAlpha = 0.2; f.EdgeColor = 'none';
% f = fill3([-50 50 50 -50], [-lim -lim -lim -lim], [5 5 15 15], [0.93 0.69 0.13 0.3]);
% f.FaceColor = 'y'; f.FaceAlpha = 0.2; f.EdgeColor = 'none';

% 3D Tube
cell_ind = 11;
timespan = 1:409;
rect = [0 0 700 850];
lim = 50;
viewangle = [389.1,330.9,211.3];
new_fig(rect, 4.5, 45, [], [], [], [], "");
interpN = 4; 
sp = 0.008;
twfactor = 10; 
zfactor = 0.5; % Factor streches (or squeezes) the z axis to prevent narrowing the tube.
tube_surf(q, cell_ind, timespan, rect, lim, viewangle, twfactor, interpN, sp, zfactor);
img = imread('..\combined_stacks_quad_Oct8_Oct23_Oct15(38)(2).tif', 150);
hold on; overlay_image(q, img, cell_ind, lim, 0, 4);
xticklabels(repmat("", [3 1]))
yticklabels(repmat("", [3 1]))
zlim([0 10*30*zfactor])
camroll(-90)
setfigprop(1300,1200,4,40)
ztickangle(-90)
set(gca, 'innerposition', [0.13,0.11,0.775,0.74], 'position', [0.13,0.11,0.775,0.74])
zticks(0:37.5:150);
zticklabels(string(0:150:600))
zl = zlabel('t (min)');
t00 = 30;
f = fill3([-lim lim lim -lim], [-lim -lim -lim -lim]+0.01, [t(t00) t(t00) t(t00+90) t(t00+90)]*60/4, ...
    [1.00,0.41,0.16]); f.FaceAlpha = 0.2; f.EdgeColor = 'none';
f = fill3([-lim -lim -lim -lim]+0.01, [-lim lim lim -lim], [t(t00) t(t00) t(t00+90) t(t00+90)]*60/4, ...
    [1.00,0.41,0.16]); f.FaceAlpha = 0.2; f.EdgeColor = 'none';
t00 = 130;
f = fill3([-lim lim lim -lim], [-lim -lim -lim -lim]+0.01, [t(t00) t(t00) t(t00+90) t(t00+90)]*60/4, ...
    [0.07,0.62,1.00]); f.FaceAlpha = 0.2; f.EdgeColor = 'none';
f = fill3([-lim -lim -lim -lim]+0.01, [-lim lim lim -lim], [t(t00) t(t00) t(t00+90) t(t00+90)]*60/4, ...
    [0.07,0.62,1.00]); f.FaceAlpha = 0.2; f.EdgeColor = 'none';

%%
% dth_smooth = smoothing(q, dth, 0.45);
dth_2_smooth = smoothing(q, dth_2, 0.85);
D_smooth = smoothing(q, D, 0.1);
th_smooth = smoothing(q, th, 0.75);

% theta vs. time (1 & 2 , 3 & 4);
init_th = initial_th(q, 1);
th_smooth2 = th_smooth + init_th;
y = squeeze(th_smooth2(:, cell_ind, :));
% Exchange 1, 4
% temp = y(1,:);  
% y(1,:) = y(4,:);
% y(4,:) = temp;
% Decrease 4rd cell by 2pi
% y(3,:) = y(3,:) + pi;
y(1,:) = y(1,:) + 2*pi;
% y(3,:) = y(3,:) + 2*pi;
y = y';
t = 0:q.T-1; t = t / 30;
% 1 & 2
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], "");
% Highlighting box
% fill([5.2 10.2 10.2 5.2], [0 0 10.3 10.3], [0.3 0.3 0.3 0.3]);
for i = 1:4
    plot(t, y(:, i), 'linewidth', 5.5, 'color', [q.colors(i, :)]);
end
xlabel("t (hour)");
ylabel("\theta (rad)")
% xlim([0 15])
xlim([0 10])
% ylim([-pi 4*pi])
% yt = [-pi 0 pi 2*pi 3*pi 4*pi]; yticks(yt);
yt = [0 2*pi 4*pi]; yticks(yt);
% ytl = {'-\pi' '0' '\pi' '2\pi' '3\pi' '4\pi'}; yticklabels(ytl)
ytl = {'0' '2\pi' '4\pi'}; yticklabels(ytl)

% theta vs. time (with different initial time)
t0 = 1;
init_th = initial_th(q, t0);
th2 = th - th(:, :, t0) + init_th;
y = squeeze(th2(:, cell_ind, t0:end));
% Exchange 1, 4
% temp = y(1,:);  
% y(1,:) = y(4,:);
% y(4,:) = temp;
y_smooth = nan(size(y));
for i = 1:4
    yy = y(i, :);
    yy = yy(~isnan(yy));
    xx = 1:length(yy);
    f = fit(xx', yy', 'smoothingspline', 'smoothingparam', 0.75);
    y_smooth(i, 1:length(yy)) = f(xx);
end
y = y';
y_smooth = y_smooth';
t = t0:q.T; t = t / 30;
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], "");
% Shading 
t00 = 30;
f = fill([t(t00) t(t00+90) t(t00+90) t(t00)], [-30 -30 30 30], [1.00,0.41,0.16]); f.FaceAlpha = 0.25; f.EdgeColor = 'none';
t00 = 130;
f = fill([t(t00) t(t00+90) t(t00+90) t(t00)], [-30 -30 30 30], [0.07,0.62,1.00]); f.FaceAlpha = 0.25; f.EdgeColor = 'none';
for i = 1:4
    plot(t, y_smooth(:, i), 'linewidth', 6, 'color', q.colors(i,:));
%     plot(t, y(:,i), '.', 'markersize', 20, 'color', q.colors(i,:)); %
%     Smoothing is not much distinguished
end
% xlabel("t (hour)");
ylabel("\theta (rad)")
% xlim([7 12])
% ylim([-pi 2*pi])
% yt = [-pi 0 pi 2*pi 3*pi 4*pi]; yticks(yt);
yt = -pi:pi:2*pi; yticks(yt)
yt = {'-\pi', '0', '\pi', '2\pi'}; yticklabels(yt)
% ylim([-pi 5.7])
ylim([-2*pi 2*pi])
xlim([0 10])
r = refline(0, 0); r.LineWidth = 2.0; r.LineStyle = '--'; r.Color = [0.3 0.3 0.3 1];
xt = 0:2.5:15; xticks(xt);
xt = xt * 60; xticklabels(xt)
setfigprop(943, 457, 4, 40)
box on

% theta vs. omega
time = (0:q.T-1)/30;
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], "");
for i = 3
% subplot(4,1,i)
yyaxis left;
plot(time, squeeze(th_smooth(i, cell_ind, :)), 'linewidth', 6, 'color', q.colors(i, :));
ylabel('\theta (rad)')
yt = 0:pi:2*pi; yticks(yt);
ytl = {'0', '\pi', '2\pi'}; yticklabels(ytl);
yyaxis right;
plot(time(1:end-dt), squeeze(dth_2_smooth(i, cell_ind, :)*30), 'linewidth', 3, 'color', [0.3 0.3 0.3 0.9]); hold on 
scatter(time(1:end-dt), squeeze(dth_2(i, cell_ind, :)*30), 50, 'markerfacecolor', [0.3 0.3 0.3], 'markeredgecolor', 'none');
r = refline(0,0); r.LineWidth = 2; r.LineStyle = '--';
ylabel('\omega (rad/h)'); 
set(gca,'fontsize', 45, 'linewidth', 4.5)
% xticks([0:3:15])
% xticklabels(["", "", "", ""])
yt = -pi:pi:3*pi; yticks(yt);
ytl = {'-\pi', '0', '\pi', '2\pi', '3\pi'}; yticklabels(ytl);
xlim([0 15])
ylim([-5 10])
box off
end
xlabel('t (hour)')

% omega vs. D
ms = 150; lw = 6;
time = (0:q.T-1)/30;
for i = 1:4
new_fig([50 50 500 500], 4.5, 45, [], [], [], [], "");
% subplot(2,2,i)
yyaxis left; 
plot(time(1:end-dt), squeeze(abs(dth_2_smooth(i, cell_ind, :))) * 30 / 60, 'linewidth', lw, 'color', [q.gra 0.9]);
scatter(time(1:end-dt), squeeze(abs(dth_2(i, cell_ind, :)*30 / 60)), ms, 'markerfacecolor', q.gra, ...
    'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', q.gra)
ylim([0-0.04 0.2+0.04])
% ylim([0 3*pi])
% ylabel('|\omega| (rad/h)')
% yt = 0:pi:3*pi; yticks(yt);
% yt = string(0:3) + "\pi"; yt(1) = '0'; yticklabels(yt);
yticks([0 0.2])
xt = 0:2.5:10; xticks(xt);
xticklabels(string(xt*60));
yyaxis right; 
plot(time, squeeze(D_smooth(i, cell_ind, :)), 'linewidth', lw, 'color', q.gre);
scatter(time, squeeze(D(i, cell_ind, :)), ms, 'markerfacecolor', q.gre, ...
    'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', q.gre)
yticks([10 30])
% ylabel('\itd (\mum)')
% xlim([5 15])
% xlim([0-0.2 10+0.2])
xlim([0 10])
ylim([4 35])
% xlim([0 15])
box on
% xlabel('t (min)')
setfigprop(918, 305, 4, 40)
% set(gca,'innerposition', [0.13,0.4,0.7453,0.6])
set(gca, 'units' , 'pixels', 'innerposition', [120.34,122.867,702.5,182.8]);
end


% Histograms
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(periodD*60, 'binwidth', 30, 'normalization', 'probability');
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0.1; h.LineWidth = 1.5;
xlabel('\tau_{\itd} (min)')
xticks(120:120:480)
xlim([20 400])
scatter(nanmean(periodD(:)*60), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 470, 4, 40)

new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(pcorrel, 'binwidth', 0.1, 'normalization', 'probability');
h.FaceColor = [0.00,0.45,0.74]; h.EdgeAlpha = 0.1; h.LineWidth = 1.5;
xlabel('\rho')
xlim([-0.6 0.6])
scatter(nanmean(pcorrel(:)), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 440, 4, 40)

% new_fig([50 50 500 600], 3.5, 35, [], [], [], [], "");
% h = histogram(mpc, 'binwidth', 0.1, 'normalization', 'probability');
% h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0;
% xlabel('c_{pair}')

%% Quadorder & transition
% Import 'quadorder(2).csv' as table
quadorder = quadorder2{:, :};
figure;
for i = 1:40
    subplot(5, 8, i)
    plot(quadorder(i, :))
    title(i)
end

%%
% States (total 6)
new_fig([50 50 500 500], 1.5, 15, [], [], [], [], "");
h = histogram(quadorder(:), 'normalization', 'probability');
h.FaceColor = [0.93,0.69,0.13]; h.EdgeAlpha = 0.1; h.LineWidth = 2.5;
xticks(1:6); xlabel("state")
setfigprop(500, 440, 4, 40)

% Types of transitions
trans = [1 2; 1 3; 1 4; 1 5; 1 6; 2 3; 2 4; 2 5; 2 6; 3 4; 3 5; 3 6; 4 5; 4 6; 5 6];
dt = 10;
quadorder_trans = nan(q.N, q.T-dt);
for n = 1:q.N-dt
    y = squeeze(quadorder(n, :));
    for t = 1:q.T-1
        qo0 = y(t);
        qo1 = y(t+dt);
        if qo0 ~= qo1
            for i = 1:15
                ismem = ismember([qo0 qo1], trans(i,:));
                if sum(ismem) == 2
                    quadorder_trans(n, t) = i;
                    continue;
                end
            end
        else
            quadorder_trans(n,t) = 0;
        end
    end
end
            
qo = quadorder_trans(:);
ind = 1:15; ind1 = [1 10 15];
ind = setdiff(ind, ind1);
for i = 1:length(qo(:))
    if ismember(qo(i), ind)
        qo(i) = 1;
    elseif ismember(qo(i), ind1)
        qo(i) = 2;
    end
end
qo(qo == 0) = [];
qo(isnan(qo)) = [];

%%
% Transition types
new_fig([50 50 500 600], 3.5, 35, [], [], [], [], "");
% h = histogram(quadorder_trans(quadorder_trans ~= 0 & ~isnan(quadorder_trans)), ...
%     'normalization', 'probability');
% h = histogram(qo, ...
%     'normalization', 'probability');
% h.FaceColor = [1.00,0.41,0.16]; h.EdgeAlpha = 0.1; h.LineWidth = 2.5;
h = histcounts(qo); h = h / sum(h);
b = bar(1:2, h);
b.BarWidth = 0.4; b.FaceAlpha = 0.5; b.EdgeColor = 'none';
b.FaceColor = [0.85,0.33,0.10];
xticks([1,2]);
xlim([0.5 2.5])
xticklabels({"", ""})
% setfigprop(500, 440, 4, 40)
setfigprop(500, 365, 4, 40)

% set(gca,'TickLength',[0.03, 0.01])

% trans_freq = sum(quadorder_trans ~= 0 & ~isnan(quadorder_trans), 2) ... 
%     ./ sum(~isnan(quadorder_trans), 2) * 30;
% 
% % Frequency?
% new_fig([50 50 600 600], 1.5, 15, [], [], [], [], "");
% histogram(trans_freq, 'binwidth', 0.35, ...
%     'normalization', 'probability'); xlabel("Mean  # of exchange per hour")





%% Explanation figures
sc = 0.15;
x = [-1 0 1 0] * sc; 
y = [0 1 0 -1] * sc;
r = [x; y];
angle = pi/8;
M = [cos(angle) -sin(angle); sin(angle) cos(angle)];
r_prime = M * r;

% scheme = [1 2 3 4; 2 1 3 4; 3 2 1 4; 1 2 4 3];
scheme = [4 1 2 3; 4 2 1 3; 4 1 2 3; 4 3 2 1];
figure; hold on 
ms = 2000;
for sp = 1:4
    subplot(1,4,sp)
    colorscheme = q.colors(scheme(sp,:), :);
    hold on 
for i = 1:4
    color = colorscheme(i, :);
    scatter(r_prime(1,i), r_prime(2,i), ms, color, 'filled')
end
lim = 0.3;
xlim([-lim lim])
ylim([-lim lim])
setfigprop(300, 300, 3, 40)
axis square
box on
set(gca,'xtick',[])
set(gca,'ytick',[])
end
setfigprop(1100, 220, 3, 40)

%% Quadruplet schematic
factor = 1;     % Scale of x, y
rfactor = 0.95;   % Scale of x-lim, y-lim
ms = 6500;       % Markersize

% order = [1 2 3 4; 
%     2 1 3 4;
%     3 2 1 4;
%     4 2 3 1];

order = [1 2 3 4; 
    2 1 3 4; 
    1 2 3 4; 
    3 2 1 4];

x = [-1 1 1 -1];
y = [-1 -1 1 1];
r = [x; y];
angle = pi/8 + pi;
% angle = 0;
M = [cos(angle) -sin(angle); sin(angle) cos(angle)];
r = M * r * factor;

for i = 1:4
    ord = order(i, :);
% order = [1 2 3 4];
% order = [2 1 3 4];
% order = [3 2 1 4];
% order = [4 2 3 1];

figure('pos', [500 500 300 300]); hold on
colors = q.colors(ord, :);
% a = 0.5;
% colors = [q.yel; ones(1,3)*a; ones(1,3)*a; ones(1,3)*a];
% colors = colors(ord, :);
for j = 1:4
sc = scatter(r(1,j), r(2,j), ms, 'filled');
if j == i || j == 1
    sc.MarkerFaceAlpha = 1;
    sc.MarkerFaceColor = colors(j, :) * 1.0;
else
%     sc.MarkerFaceAlpha = 0.4;
%     sc.MarkerFaceColor = colors(j, :) * 0.9;
    sc.MarkerFaceAlpha = 1;
    sc.MarkerFaceColor = colors(j, :);
    
end
end
xlim(rfactor*[-2 2])
ylim(rfactor*[-2 2])
axis square
axis off
end

%% Save all figures
figHandles = findall(0,'Type','figure'); 

savefig(figHandles, "figures\figs_pair_quad_dec_21.fig");


%% Checking
y = q.pos(:, cell_ind, :, :);
y = squeeze(y);
% y(:,2,:) = - y(:,2,:);
my = squeeze(nanmean(y, 1));
t = 200:18:272;
figure;
for i = 1:5
    subplot(3,2,i)
    scatter(y(:, 1, t(i)), y(:, 2, t(i))); hold on 
    scatter(my(1, t(i)), my(2, t(i)));
    xlim([220 245])
    ylim([-190 -165])
    axis square
end

y1 = y(1,:);
y2 = y(2,:);
y3 = y(3,:);
y4 = y(4,:);
for t = 1:length(y1)
    if abs(y4(t) - y3(t)) > 2*pi
        y3(t) = y3(t) + 2*pi;
    end
end

figure; hold on 
t = (0:length(y1)-1) / 30;
plot(t, y1, 'linewidth', 6)
plot(t, y2, 'linewidth', 6)
plot(t, y3, 'linewidth', 6)
plot(t, y4, 'linewidth', 6)


%% Draw trajectory on the image 
% Read images
imgs = []; t = 1;
% for i = 0:2:240-1
for i = 0:298
    name = ['quad_154_452_crop(idx11)(2)\img_' num2str(i,'%03.f') '.png'];
    A = imread(name);
    imgs(:, :, t) = A;
    t = t + 1;
end

%% Smooth the curves
% Check cell_ind (currently 11)
quad = q.pos(1:4, cell_ind, :, 5:303);
quad = permute(quad, [4,3,2,1]);
quad = squeeze(quad);
offset = nanmean(quad, 3);
quad = quad - offset(1, :);
quad = quad + 16;
quad(:, 1, :) = quad(:, 1, :)  + 0.8;
quad(:, 2, :) = quad(:, 2, :) + 2.5;
qcm = nanmean(quad , 3);

[T, ~, N] = size(quad);
q_smooth = nan(size(quad));
for i = 1:4
    for j = 1:2
        y = quad(:, j, i); ind = find(~isnan(y));
        x = 1:length(y);
        y = y(ind); x = x(ind);
        f = fit(x', y, 'smoothingspline', 'smoothingparam', 0.01);
        q_smooth(ind, j, i) = f(x);
    end
end

%%
% View & Save overlay images
lw = 6;
ms = 50;
fig = figure('pos', [100 100 700 700]);
% for t0 = 1:120
for t0 = 100
    t = t0;
    imshow(imgs(:, :, t) / 255, 'InitialMagnification',1000); hold on;
    celln = 1;
    plot(q_smooth(t0:t, 1, celln), q_smooth(t0:t, 2, celln), 'linewidth', lw, 'color', q.colors(celln, :))
    plot(q_smooth(t, 1, celln), q_smooth(t, 2, celln), '.', 'markersize', ms, 'color', q.colors(celln, :))
    celln = 2;
    plot(q_smooth(t0:t, 1, celln), q_smooth(t0:t, 2, celln), 'linewidth', lw, 'color', q.colors(celln, :))
    plot(q_smooth(t, 1, celln), q_smooth(t, 2, celln), '.', 'markersize', ms, 'color', q.colors(celln, :))
    celln = 3;
    plot(q_smooth(t0:t, 1, celln), q_smooth(t0:t, 2, celln), 'linewidth', lw, 'color', q.colors(celln, :))
    plot(q_smooth(t, 1, celln), q_smooth(t, 2, celln), '.', 'markersize', ms, 'color', q.colors(celln, :))
    celln = 4;
    plot(q_smooth(t0:t, 1, celln), q_smooth(t0:t, 2, celln), 'linewidth', lw, 'color', q.colors(celln, :))
    plot(q_smooth(t, 1, celln), q_smooth(t, 2, celln), '.', 'markersize', ms, 'color', q.colors(celln, :))
    
    plot(qcm(t, 1), qcm(t, 2), 'x', 'markersize', 15, 'color', 'w', 'linewidth', 3)
    hold off
%     saveas(gcf, ['img_' num2str(t0) '.png']);
end

%% Overlay with trajectory
lw = 11;
ms = 80;
alpha = 0.7;
alpha1 = 1.0;
% for dt = [0 27 54 81]
for dt = [0 30 60 90]
figure;
% n = n + 1;
% sub = subplot(2, 4, n);
% ps = get(sub, 'Position');
% w = ps(3); h = ps(4);
% %%%%% NOTE: This is 4 mins interval %%%%%%%
% t00 : index of frame from the 4 mins interval sequence
% t0 : index of frame from the original sqeuence

% Cells 1, 3
% t0 = 255;  % 8.5 (hour)
% t0 = 280;
% Cells 1, 2
% t0 = 210;       % 7 (hour)
% Cells 3, 4
% t0 = 180;   % 6 (hour)
% t0 = 255;
% Cells 2, 3
% t0 = 240;

t0 = 186; % 6.2 hr
% t0 = 210; % 7 hr
% t0 = 240; % 8 hr
% t0 = 270; % 9 hr

t00 = t0 - 120 ;t00 = t00 / 2; t00 = floor(t00);

t = t00 + 0;
t = t00 + 15; % 1hr
t = t00 + 30; % 2hr
t = t00 + 45; % 3hr 
% t = t00 + 60; % 4hr

% t = t00 + 40;      % 2.7 hours (Cell 1, 3)
% t = t00 + 45; % 3 hours (Cell 3, 4)
% t = t00 + 37; % 2.5 hours (Cell 2, 3 and 1, 2)

%%%%% Starting time %%%%
% t0 = 189; % 6hour
t0 = 280; % 10hour
t00 = t0 - 154 + 1;
% dt = 80;
t = t00 + dt; % 3hours


% t00 = 30; % 1 h
% t = t00 + dt;

t00 = 130; % 5 h
t = t00 + dt;

%%%%                %%%%

cq = nanmean(q_smooth, 3);
% t0 = 55; t = 120;
% t0 = 120; t = 120;
imshow(imgs(:, :, t) / 255, 'InitialMagnification',1000); hold on;
for celln = 4:-1:1
% plot(q_smooth(t0, 1, celln), q_smooth(t0, 2, celln), '.', 'markersize', ms, 'color', [q.colors(celln, :) alpha])
plot(q_smooth(t00:t, 1, celln), q_smooth(t00:t, 2, celln), 'linewidth', lw, 'color', [q.colors(celln, :)*0.9 alpha])
plot(q_smooth(t, 1, celln), q_smooth(t, 2, celln), '.', 'markersize', ms, 'color', [q.colors(celln, :) alpha1])
end
% plot(cq(t00:t, 1), cq(t00:t, 2), 'linewidth', 4, 'color', [0.9 0.9 0.9 0.6])
plot(cq(t, 1), cq(t, 2), 'x', 'linewidth', 4, 'markersize', 15, 'color', [0.9 0.9 0.9 alpha])

end

%%
t = 120;
y = quad;
y(:, 2, :) = y(:, 2, :);
y(:, 1, :) = y(:, 1, :);
figure; 
imshow(imgs(:, :, t) / 255);
hold on
yy = squeeze(y(t,:,:));
scatter(yy(1,1), yy(2,1), 'markerfacecolor', 'r')
scatter(yy(1,2), yy(2,2), 'markerfacecolor', 'b')
scatter(yy(1,3), yy(2,3), 'markerfacecolor', 'y')
scatter(yy(1,4), yy(2,4), 'markerfacecolor', 'g')










%% Checking the index of a quad (overlay with texts)
figure; 
t = 154;
img = imread('..\combined_stacks_quad_Oct8_Oct23_Oct15(38).tif', t);
imshow(img);
x = squeeze(pos_quad(:, :, 1, t));
y = squeeze(pos_quad(:, :, 2, t));
hold on 
for i = 1:40
    for j = 1:4
        xx = x(j, i);
        yy = y(j, i);
        scatter(xx, yy, 10);
        t = text(xx, yy, num2str(i));
        t.Color = q.colors(j, :);
    end
end

%%








































%% 0. Calculate the data
position = pos_quad;
T = size(position, 4);

mp = nanmean(position, 1);
dv = position - mp;
adv = absvec(dv, 3);
adv(adv==0) = 1j;
du = dv./adv;
du = real(du);

% Angular velocity
dth = nan(size(du));
dth = squeeze(dth(:,:,1,1:end-1));


for t=1:T-1
    u1 = du(:,:,:,t);
    u2 = du(:,:,:,t+1);
    crossu = u1(:,:,1).*u2(:,:,2) - u1(:,:,2).*u2(:,:,1);
    crossu = - crossu;
    ang = asin(crossu);
    dth(:,:,t) = ang;
%     waitbar(t/T)
end


sc = 1000 / 340.5;
mpos = nanmean(position, 1);
D = position - mpos;
D = squeeze(absvec(D, 3)) * sc;



x_unit = zeros(size(du));
x_unit(:,:,1,:) = 1;
% cross = crossp(du, x_unit, 3);
cross = du(:, :, 1, :).*x_unit(:, :, 2, :) - du(:, :, 2, :).*x_unit(:, :, 1, :);
dot = du(:, :, 1, :).*x_unit(:, :, 1, :) + du(:, :, 2, :).*x_unit(:, :, 2, :);
ang = sign(cross).*acos(dot);
ang = squeeze(ang);
dth = cat(3, ang(:,:,1), dth);

% dth = cat(3, dth, zeros(5, 40, 1));

th = nan(size(dth));
for i=1:size(dth,2)
    for j=1:size(dth,1)
        a = dth(j,i,:);
        a = squeeze(a);
        a = a(~isnan(a));
        ac = cumsum(a);
        th(j,i,1:length(ac)) = ac; % Start time doesn't matter much
    end
end

%% 1. FFT
%% Subtract heavily smoothed D
newD = D;
for i=1:size(D,2)
    for j=1:size(D,1)
        y = squeeze(D(j,i,:));
        IND = find(~isnan(y));
        x = 1:length(y); x = x';
        x = x(IND);
        y = y(IND);
        if ~isempty(y)
            f = fit(x, y, 'smoothingspline', 'SmoothingParam', 3e-6);
            newD(j, i, IND) = newD(j, i, IND) - reshape(f(x), 1, 1, length(IND));
        end
    end
end

newD = reshape(newD, [], size(newD, 3));

% Check 
figure; 
for i=1:size(newD, 1)
    subplot(10, 10, i)
    plot(squeeze(newD(i,:)))
end
%% Subtract heavily smoothed dth
newDth = dth;
for i=1:size(dth, 2)
    for j=1:size(dth, 1)
        y = squeeze(dth(j, i,:)); 
        IND = find(~isnan(y));
        x = 1:length(y); x = x';
        x = x(IND);
        y = y(IND);
        if ~isempty(y)
            f = fit(x, y, 'smoothingspline', 'SmoothingParam', 3e-6);
            newDth(j, i, IND) = newDth(j, i, IND) - reshape(f(x), 1, 1, length(IND));
        end
    end
end

newDth = reshape(newDth, [], size(newDth, 3));

% Check 
figure; 
for i=1:size(newDth, 1)
    subplot(10, 10, i)
    plot(squeeze(newDth(i,:)))
end



%% FFT subplots
figure;
N = size(newD, 1);
T = 2/60; % Interval (hours)
Fs = 1/T;
maxT = nan(N ,1);
j = 1;
for i=1:1:N
    
%     dd = squeeze(D(i, :));
    dd = squeeze(newD(i, :));
    dd = dd(~isnan(dd));
    if ~isempty(dd)
        dd = dd(1:2*floor(length(dd)/2));
        L = length(dd);
%         t = (0:L-1)*T;
        Y = fft(dd);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1); 
        f = Fs*(0:(L/2))/L;
        t = 1./f;
        
        % Find peak
%         diffP1 = diff(P1);
%         diffP1 = sign(diffP1);
%         diffP12 = diff(diffP1);
%         peakInd = find(diffP12 == -2);


        
        % Smooth the curve and get the peak
        f = fit(t(2:end)', P1(2:end)', 'smoothingspline', 'SmoothingParam', 0.90);
        f = f(t(2:end));
        ind = find(max(f) == f) + 1;
        maxT(i) = t(ind);
        


%         subplot(10, 6, j)
%         yyaxis left
%         plot(t, P1, '.-', 'markersize', 10, 'linewidth', 2)     
%         
%         hold on 
%         plot(t(2:end), f, 'g')
%         plot(t(ind), 0, 'kx', 'linewidth', 2, 'markersize', 8)
%         
%         yyaxis right
%         plot((0:L-1)*T, dd)
%         ylabel('r_g')
        title(['T_r_g= ' num2str(round(t(ind),3))])

   
        xlim([0 15])
        
        j = j + 1;
    end
    
end

edges = 0:0.5:7;
H = zeros(1, length(edges)-1);
for i = 1:100
    a = maxT(~isnan(maxT));
    r = randperm(length(a), 57);
    maxT_randperm = maxT(r);
    h = histcounts(maxT_randperm, edges);
    H = H + h;
end
H = H / 100;

%% Display the oscillation period histogram
figure('pos', rect); h = bar(edges(1:end-1)+0.25, H);
h.BarWidth = 1; h.FaceColor = 'k'; h.FaceAlpha = 0.7; h.EdgeAlpha = 0;
% figure; histogram(maxT_randperm, 'binwidth', 0.5)
xlabel('T (h)')
% title('Oscillation period of r_g')
set(gca, 'fontsize', 35, 'linewidth', 3.5)
xlim([0 7])
ylim([0 12])
box off


%% 2. The angles vs time

% q = [11, 30, 25];
% q = 11;

% Set the initial time
t0 = 1;
dth_11 = dth(:, 11, t0:end);
x_unit = zeros(size(du));
x_unit(:,:,1,:) = 1;
% cross = crossp(du, x_unit, 3);
cross = du(:, :, 1, :).*x_unit(:, :, 2, :) - du(:, :, 2, :).*x_unit(:, :, 1, :);
dot = du(:, :, 1, :).*x_unit(:, :, 1, :) + du(:, :, 2, :).*x_unit(:, :, 2, :);
ang = sign(cross).*acos(dot);
ang = squeeze(ang);
A = ang(:, :, t0);
dth_11 = dth_11(:, :, t0:end);
dth_11 = cat(3, A, dth_11);

sz = size(dth_11);
th_11 = nan(sz);
for i=1:size(position,2)
    for j=1:size(position,1)
        a = dth_11(j,i,:);
        a = squeeze(a);
        a = a(~isnan(a));
        ac = cumsum(a);
        th_11(j,i,1:length(ac)) = ac; % Start time doesn't matter much
    end
end


% Every quad
T = size(th, 3);
figure; hold on
ii=  1;
for qi = 1:40
% for i=1:10
    subplot(5, 8, ii)
    hold on 
%     plot((0:T-1)*2/60, squeeze(nanmean(th(:, qi, :))), '-', 'linewidth', 6, 'color', [0.3 0 1 0.5])
    for j = 1:size(th,1)
    %         plot((0:T-1)*2/60, squeeze(th(j,pi,:)), 'linewidth', 4, 'color', 'k'); hold on 
        plot((0:T-1)*2/60, squeeze(th(j,qi,:)), 'linewidth', 2); hold on 
%             ylim([-30 30])
        ylim([-15 15])
        title(num2str(qi))
        xlim([0 14])
%         xticklabels(string(0:2:12))
%         r = refline(0, 0); r.LineStyle = ':'; r.Color = [0 0 0];
%         set(gca, 'fontsize', 25, 'linewidth', 3)
%         box off
%         ylabel('\theta (rad)')
    end
    ii = ii + 1;   
end

%%
% Selected samples
q = [11 30 25];
T = size(th, 3);
figure; hold on
ii=  1;
for qi = 30
% for i=1:10
%     subplot(3, 1, ii)
    fill([0 5 5 0], [-30 -30 30 30], 'b')
    fill([5 10 10 5], [-30 -30 30 30], 'r')
    for j = 1:4
        t = 1:558;
        y = squeeze(th(j, qi,t));
        
        t = t(~isnan(y));
        y = y(t);
        
    
        f = fit(t', y, 'smoothingspline', 'smoothingparam', 0.001);
        
        y = f(t);
        
        plot(t*2/60, y, 'linewidth', 6); hold on 
        
%         plot((0:T-1)*2/60+2, squeeze(th(j,qi,:)), 'linewidth', 6); hold on 
%             ylim([-30 30])
%         ylim([-15 15])
%         title(num2str(qi))
%         xlim([2 16])
%         xticklabels(string(0:2:12))
        r = refline(0, 0); r.LineStyle = ':'; r.Color = [0 0 0];
%         set(gca, 'fontsize', 25, 'linewidth', 3)
        box off
        ylabel('\theta (rad)')
    end
    ii = ii + 1;   
end

ylim([-pi 4*pi])

xt = t(1):120:t(end);
xtl = string((0:120:800)/30);
% xticks(xt);
% xticklabels(xtl);

yt = [-6*pi : 2*pi : 6*pi]; 
yticks(yt);
ytl = {'-6\pi' '-4\pi' '-2\pi' '0' '2\pi' '4\pi' '6\pi'}; 
yticklabels(ytl)

set(gca, 'fontsize', 35, 'linewidth', 3.5) 
ylabel('\theta')
xlabel('hour')

%%
% N = {};
% N{1} = [1 27 111 195 247 260 305 355];
% N{2} = [1 39 128 200];
% N{3} = [1 12 23 47 102 220 276 324];
% N{4} = [1 139 150 185 221 238 252 266];
% N{5} = [1 111 125 159 281];
% N{6} = [1 54 167 258 526];
% N{7} = [1 74 137 160 218 303 357 469];
% N{8} = [1 32 140 247 316 351];
% N{9} = [1 21 80 131 158 181 232];
% N{10} = [1 16 107 174 223];
% N{11} = [1 33 79 94 174 250 333 452];
% N{12} = [1 19 25 113 151 178];
% N{13} = [1 80 151 281];
% N{14} = [1 34 115 142 184];
% N{15} = [1 54 154 223 237 385];
% N{16} = [1 205 440];
% N{17} = [1 31 108 144 191 331];
% N{18} = [1 30 87 103 190 238 281];
% N{19} = [1 118 139 230 316];
% N{20} = [1 24 107 184 237 359 395];
% N{21} = [1 34 66 97 172 192 272 318];
% N{22} = [1 43 251 379];
% N{23} = [1 59 115 145 229 287];
% N{24} = [1 30 109 162 187 255 300 337 440];
% N{25} = [1 73 138 169 226 270 294 340 398 424 499];
% N{26} = [1 105 120 160 223 262 375];
% N{27} = [1 28 82 150 238 338 391 430];
% N{28} = [1 27 45 64 92 141 158 213 273 309 367 417];
% N{29} = [1 31 36 57 59 73 96 106 121 142 164 172 183 223 259 296 356 367 392 402 418 459 475 483 501];
% N{30} = [1 141 150 196 216 247 314 344  433];
% N{31} = [1 9 14 61 93 109 118 139 147 169 176 203 295 347];
% N{32} = [23 67 92 122 146 189 290 416];
% N{33} = [1 27 81 126];
% N{34} = [1 53 106 138 199 272 439];
% N{35} = [1 15 23 59 130 149 167];
% N{36} = [1 60 112 193 263 327 370 408];
% N{37} = [1 48 67 105 124 132 150 161 166 180 213 221 236 258];
% N{38} = [1 18 27 57 94 127 353];
% N{39} = [1 18 79 95 123 149 177 194 203 216 220 236 242 249 256 282 285 301 324 349 373];
% N{40} = [1 70 152 240 319 330 391 415 488];

N = {};
N{1} = [1 27 111 195 247 260 305 355];
N{2} = [1 39 128 200];
N{3} = [1 12 23 47 102 220 276 324];
N{4} = [1 139 150 185 221 238 252 266];
N{5} = [1 37 47 63 111 125 141 154 196 281];
N{6} = [1 54 167 183 258 526];
N{7} = [1 74 137 160 218 303 357 372 469];
N{8} = [1 32 140 247 316 351];
N{9} = [1 21 80 131 158 181 232];
N{10} = [1 16 107 174 223];
N{11} = [1 33 79 94 174 250 333 452];
N{12} = [1 18 22 111 137 151 178];
N{13} = [1 80 151 240 281];
N{14} = [1 34 115 142 184];
N{15} = [1 54 154 223 237 385];
N{16} = [1 232 272 440];
N{17} = [1 31 108 144 191 331];
N{18} = [1 30 87 103 162 187 238 281];
N{19} = [1 118 139 230 316];
N{20} = [1 24 107 184 237 359 395];
N{21} = [1 34 66 97 172 192 272 318];
N{22} = [1 43 251 379];
N{23} = [1 59 115 145 229 287];
N{24} = [1 30 109 162 187 255 300 337 389 440];
N{25} = [1 73 138 169 226 270 294 340 398 424 499];
N{26} = [1 105 120 160 223 262 375];
N{27} = [1 28 82 150 238 338 391 430];
N{28} = [1 45 71 92 141 158 213 273 309 367 417];
N{29} = [1 59 73 96 121 164 206 259 296 356 369 393 418 487 515];
N{30} = [1 141 216 247 433];
N{31} = [1 61 110 169 259 347];
N{32} = [23 67 92 122 189 416];
N{33} = [1 27 81 126];
N{34} = [1 53 106 138 272 439];
N{35} = [1 15 23 59 130 149 167];
N{36} = [1 112 193 263 327 408];
N{37} = [1 48 67 105 150 213 236 258];
N{38} = [1 18 27 94 353];
N{39} = [1 18 79 95 123 177 280 293 324 350 373];
N{40} = [1 70 152 240 295 319 330 391 415 488];



q = CellCluster(pos_quad);
dth = ang_vel(q);
D = cell_dist(q);
th = ang(q);

th(:, 33, 130:end) = nan;
th(:, 40, 319:end) = nan;

% Average the angle (th)
% mth = squeeze(mean(th, 1));
mth = nan(size(th, 2), size(th, 3));
for n = 1:40
    a = th(:, n, :);
    a = squeeze(a);
    c = [];
    for i = 1:5
        b = a(i, :);
        if sum(isnan(b)) ~= length(b)
            c = [c; b];
        end
    end
    d = nanmean(c, 1);
    mth(n, 1:length(d)) = d;
end



% q = [11 30 25];
T = size(th, 3);
figure; hold on
ii=  1;
for qi = 1:40
    subplot(5, 8, ii); hold on 
    shade = shadedErrorBar((0:T-1)*2, mth(qi, :), nanstd(mth, 0, 1));
    shade.edge(1).LineStyle = 'none'; shade.edge(2).LineStyle = 'none';
    shade.patch.FaceColor = q.blu;
    
    plot((0:T-1)*2, mth(qi, :), 'linewidth', 2)
    for n = N{qi}
        plot(n*2, mth(qi,n), 'x', 'linewidth', 2, 'markersize', 6)
    end
%     title(num2str(qi))
    r = refline(0, 0); r.LineWidth = 1; r.LineStyle = '--';
    ylabel('\theta')
    xlabel('min')
    ylim([-20 20])
    box on
    ii = ii + 1;   
end

%%
slopeq = {};
for c = 1:40
    dy = diff(mth(c, N{c})); % In radians
    dx = diff(N{c})/30; % in hours
    slopeq{c} = dy./dx;
end

% Mean of each slopes
meanSlopeq = [];
for c = 1:40
    s = slopeq{c};
    meanSlopeq(c) = nanmean(abs(s));
end

%% Histogram of angular speed
% figure('pos', rect); 
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(meanSlopeq, 'binwidth', 0.2, 'normalization', 'probability');
% h.FaceColor = 'k'; h.FaceAlpha = 0.7; h.EdgeAlpha = 0;
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0;
xlabel('\omega (rad/h)')
% ylim([0 10])
xlim([0 2])
% title('angular speed')
% set(gca, 'fontsize', 35, 'linewidth', 3.5)
box off

%% Display the rotation period (inverse of angular velocity)
% rect = [100 100 600 600];
% figure('pos', rect); 
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(2*pi ./ meanSlopeq * 60, 'BinWidth', 140, 'normalization', 'probability'); 
% h.FaceColor = 'k'; h.FaceAlpha = 0.7; h.EdgeAlpha = 0;
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0.1; h.LineWidth = 1.5;
xlabel('{\itT} (min)')
% xlim([50 1500])
xlim([50 1500])
xticks(240:480:3000)
% xticks(240:720:3000)
setfigprop(500, 440, 4, 40)
% ylabel('n')
% xlim([2 20])
% ylim([0 20])
% title('angular speed')
% set(gca, 'fontsize', 35, 'linewidth', 3.5)
scatter(nanmean(2*pi ./ meanSlopeq * 60), 0.0, 350, [1 0 0], '^', 'filled')
box off

%% Now, I can 'flatten' the theta vs. time
% Make all negative slopes from 'Mnew' positive, and cumulative sum to get
% the angles.

mdth = diff(mth, 1, 2); 
thFlat = nan(size(mth));
time = nan(size(mth));
for c=1:40
%     s = diff(mth(c,N{c}));
%     ssign = sign(s(1));
%     as = ssign * abs(s);
%     theta = cumsum(as);
%     theta = [0 theta];
%     thFlat(c, 1:length(theta)) = theta;
%     
%     t = (N{c}-1) / 30; % hour
%     time(c, 1:length(t)) = t;
    s = diff(mth(c,N{c}));
    dtheta = [];
    ssign = sign(s(1));
    for si = 1:length(s)
        if sign(s(si)) ~= ssign 
            w = - mdth(c, N{c}(si) : N{c}(si+1));
        else
            w = mdth(c, N{c}(si) : N{c}(si+1));
        end
        dtheta = [dtheta w];
    end
    w = w / 30;
%     as = ssign * abs(s);
    theta = cumsum(dtheta);
    theta = [0 theta];
    thFlat(c, 1:length(theta)) = theta;
    t = (N{c}-1) / 30; % hour
    time(c, 1:length(t)) = t;
end


% Display the result
T = size(mth, 2);
figure; hold on 
for c=1:40
    subplot(7,9,c)
    plot(time(c,:), thFlat(c,:), 's-', 'linewidth', 1.5); hold on
%     plot((0:T-1)/30, mth(c,:), 'linewidth', 1.5)
    plot((0:T-1)/30, thFlat(c, :))
    for m = N{c}
        plot((m-1)/30, mth(c,m), 'x', 'linewidth', 2)
    end
    ylim([-40 40])
%     title(num2str(c))
end

%% All new angles on a graph
figure('pos', [0,0,693,479]); hold on 
for c=1:40
    plot((0:558-1)/30, thFlat(c, :), 'linewidth', 2.4, 'color', [rand(1,3) 0.6]);
end
xlim([0 20])
r = refline(0, 0); r.Color = 'k'; r.LineWidth = 3; r.LineStyle = '--';
xlabel('t (min)')
ylabel('<\theta> (rad)')
ylim([-60 60])
set(gca, 'fontsize', 37, 'linewidth', 3.7)
yt = [-16*pi : 4*pi : 16*pi]; yticks(yt);
ytl = {'-16\pi' '-12\pi' '-8\pi' '-4\pi' '0' '4\pi' '8\pi' '12\pi' '16\pi'}; yticklabels(ytl)
xt = 0:5:20; xticks(xt)
xt = xt * 60; xticklabels(xt)
ylim([-8*pi 8*pi])
setfigprop(1200, 800, 4, 40)
box off

%% 3. Display the dth(left), D(right) vs time : QUADS
T = size(dth, 3);
i = 1;
len = 20;
qq = [1, 3, 4]; % Cell index 
q = [11, 30, 25];
% col = {'r', 'b', [0 0.7 0], [0.93 0.69 0.13]};
figure; i = 2;
for qi = 11
    for i = 1:4
        subplot(4,1,i)
        %     f = abs(dth(pi, :)) .* D(pi,1:end-1) * sc;  % Tangential speed
        f = squeeze(abs(dth(i, qi, :))); % Only the 'w'
        %     plot(0:T-2, squeeze(abs(dth(p, n, :)))); ylabel('|\omega|D')
        f_new = nan(T-len,1);
        for t=1:T-len-1
            f_new(t) = mean(f(t:t+len));
        end

        yyaxis left
        %     plot((0:T-len-1)/(T-len-1)*(T-1)*2/60, squeeze(f_new) * 30, 'linewidth', 4); ylabel('v_{\theta}(\mum/h)')
        plot((0:T-len-1)/(T-len-1)*(T-1)*2/60, squeeze(f_new) * 30, 'linewidth', 4); ylabel('\omega (rad/h)')
        ylim([0 6])
        yticks(0:3:6)

        yyaxis right
        plot((0:T-1)*2/60, squeeze(D(i, qi ,:)) , 'linewidth', 4); ylabel('r_g(\mum)')
        ylim([0 50])
        xlim([0 15])
        %     plot((0:T)*2/60, shape, 'linewidth', 4); ylabel('shape index')
        yticks(0:20:40)
        %     xticklabels(string(0:2:12))
        i = i + 1;
        set(gca, 'fontsize', 20, 'linewidth', 2)
        box off
    end
end
xlabel('hour')


%% 4. Pearson correlation (histogram)
% X = D(:, 1:end-1);
% Y = abs(dth(:, :)) * 30; % Rad per hour

X = newD;
Y = abs(newDth) * 30; % Rad per hour

mX = nanmean(X, 2);
mY = nanmean(Y, 2);
A = (X - mX) .* (Y - mY);
A = nansum(A, 2);
B = sqrt(nansum((X - mX).^2, 2) .* nansum((Y - mY).^2, 2));
r = A ./ B;

% Verify the graphs 
figure; 
for i=1:57
    subplot(10, 6, i)
    scatter(X(i,:), Y(i,:), 2)
    title(num2str(r(i)))
end


% Randomly permuate 57 cells out of 200 cells 100 times, average.
edges = [-0.5250 -0.4500 -0.3750 -0.3000 -0.2250 -0.1500 -0.0750 0 0.0750 0.1500];
H = zeros(1, length(edges)-1);
for i = 1:100
    a = maxT(~isnan(r));
    rp = randperm(length(a), 57);
    r_randperm = r(rp);
    h = histcounts(r_randperm, edges);
    H = H + h;
end
H = H / 100;


%% Draw histogram
figure('pos', rect); h = bar(edges(1:end-1)+0.075/2, H);
h.BarWidth = 1; h.BarWidth = 1; h.FaceColor = 'k'; h.FaceAlpha = 0.7; h.EdgeAlpha = 0;
xlim([-0.6 0.2])
ylim([0 12])
set(gca, 'fontsize', 35, 'linewidth', 3.5)
xlabel('\rho')
% title('Pearson correlation coefficient')
box off


%% Backup
N = {};
N{1} = [1 27 111 195 247 260 305 355];
N{2} = [1 39 128 200];
N{3} = [1 12 23 47 102 220 276 324];
N{4} = [1 139 150 185 221 238 252 266];
N{5} = [1 111 125 159 281];
N{6} = [1 54 167 258 526];
N{7} = [1 74 137 160 218 303 357 469];
N{8} = [1 32 140 247 316 351];
N{9} = [1 21 80 131 158 181 232];
N{10} = [1 16 107 174 223];
N{11} = [1 33 79 94 174 250 333 452];
N{12} = [1 19 25 113 151 178];
N{13} = [1 80 151 281];
N{14} = [1 34 115 142 184];
N{15} = [1 54 154 223 237 385];
N{16} = [1 205 440];
N{17} = [1 31 108 144 191 331];
N{18} = [1 30 87 103 190 238 281];
N{19} = [1 118 139 230 316];
N{20} = [1 24 107 184 237 359 395];
N{21} = [1 34 66 97 172 192 272 318];
N{22} = [1 43 251 379];
N{23} = [1 59 115 145 229 287];
N{24} = [1 30 109 162 187 255 300 337 440];
N{25} = [1 73 138 169 226 270 294 340 398 424 499];
N{26} = [1 105 120 160 223 262 375];
N{27} = [1 28 82 150 238 338 391 430];
N{28} = [1 45 71 92 141 158 213 273 309 367 417];
N{29} = [1 59 73 96 121 164 259 296 356 376 418 487 515];
N{30} = [1 141 216 247 433];
N{31} = [1 110 169 347];
N{32} = [23 67 92 122 189 416];
N{33} = [1 27 81 126];
N{34} = [1 53 106 138 272 439];
N{35} = [1 15 23 59 130 149 167];
N{36} = [1 112 193 263 327 408];
N{37} = [1 48 67 105 150 213 236 258];
N{38} = [1 18 27 94 353];
N{39} = [1 18 79 95 123 177 285 324 350 373];
N{40} = [1 70 152 240 319 330 391 415 488];




































