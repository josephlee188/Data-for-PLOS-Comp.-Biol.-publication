trange = 1000:4000;
pos_quad_simul = permute(pos(trange, :, :, :), [3 4 2 1]);
pos_quad_simul = cat(1, pos_quad_simul, nan(1, 200, 2, length(trange)));
qs = CellCluster(pos_quad_simul);
% qs.scale = 1;
qs.scale = 3;

dth = ang_vel(qs);
D = cell_dist(qs);
th = ang(qs);
qc = quad_correlation(qs);
dt = 5;
dth_2 = nan(5, qs.N, qs.T-dt);
for t = 1:qs.T-dt
    dth_2(:, :, t) = th(:, :, t+dt) - th(:, :, t);
end
dth_2 = dth_2 / dt;

newD = subtract_smoothing(qs, D, 3e-6);
newDth = subtract_smoothing(qs, dth, 3e-6) * 30; % rad per hour

periodD = get_period_fft(qs, newD); % In hour
periodDth = get_period_fft(qs, newDth); % hour

pcorrel = pearson_correlation(qs, newD, abs(newDth));


%% Quads
dth_2_smooth = smoothing(qs, dth_2, 0.85);
D_smooth = smoothing(qs, D, 0.1);
th_smooth = smoothing(qs, th, 0.75);

% theta vs. time
cell_ind = 2;
% t0 = 2000;
% t0 = 2000; 
t0 = 1; t0 = 2400;
init_th = initial_th(qs, t0);
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
y_smooth(:, 3) = y_smooth(:, 3) + 2*pi;
t = t0:qs.T; 
t = t / 30;
% 1 & 2
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], "");
% Shading background
f = fill([0 3 3 0], [-20 -20 20 20], [0.07,0.62,1.00]); f.FaceAlpha = 0.35; f.EdgeColor = 'none';
f = fill([3.66 6.66 6.66 3.66], [-20 -20 20 20], [1.00,0.41,0.16]); f.FaceAlpha = 0.35; f.EdgeColor = 'none';
for i = 1:4
    plot(t - t(1), y_smooth(:, i), 'linewidth', 6, 'color', qs.colors(i,:));
end
% xlabel("t");
ylabel("\theta (rad)")
% ylim([-3*pi 3*pi])
% yt = [-6*pi:3*pi:6*pi]; yticks(yt);
ylim([-17.62 3.6])
yt = [-4*pi:2*pi:4*pi]; yticks(yt);
ytl = string(-4:2:4) + "\pi"; ytl(3) = "0";
yticklabels(ytl);
r = refline(0, 0); r.LineWidth = 3; r.LineStyle = '--'; r.Color = [0.3 0.3 0.3 1];
xlim([0 10])
xt = 0:2.5:15; xticks(xt);
xt = xt * 60; xticklabels(xt)
setfigprop(1000, 450, 4, 40)
% setfigprop(943, 457, 4, 40)
% set(gca, 'units', 'pixels', 'innerposition', [169,80,684,342])
box on

% omega vs. D
ms = 150; lw = 6;
timespan = t0 : (t0 + dur);
tspan = timespan - timespan(1); tspan = tspan * 2;
t = time - time(1); t = t / 30;
for i = 1:4
new_fig([50 50 400 400], 4.5, 45, [], [], [], [], "");
% subplot(4,1,i)
yyaxis left; 
xlim([0 600])
r = refline(0, 0); r.LineWidth = 2; r.LineStyle = '--'; r.Color = 'k';
plot(tspan, squeeze(abs(dth_2_smooth(i, cell_ind, timespan))) * 30 / 60, 'linewidth', lw, 'color', [qs.gra 0.9]);
scatter(tspan, squeeze(abs(dth_2(i, cell_ind, timespan)*30/60)), ms, 'markerfacecolor', qs.gra, ...
    'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', qs.gra)
yt = [0 0.2]; yticks(yt)
ylim([-0.05 0.25])
% ylim([0 4*pi])
% ylabel('|\omega| (rad/h)')
% yt = [pi 3*pi]; yticks(yt);
% ytl = {'-\pi', '0', '\pi', '2\pi', '3\pi'}; yticklabels(ytl);
% ytl = {'', '3'}; ytl = ytl + "\pi"; yticklabels(ytl);
yyaxis right; 
plot(tspan, squeeze(D_smooth(i, cell_ind, timespan) * 3), 'linewidth', lw, 'color', qs.gre); % Multiply 3 to match um scale
scatter(tspan, squeeze(D(i, cell_ind, timespan) * 3), ms, 'markerfacecolor', qs.gre, ...
    'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', qs.gre)
% ylabel('\itd (\mum)')
% xlabel('t (hour)')
% yt = 0:4:8; yticks(yt)
yt = [5 25]; yticks(yt)
ylim([0 30])
box on
% xlabel('t (min)')
% xlim([0 600])
xt = 0:150:600; xticks(xt);
setfigprop(1000, 300, 4, 40)
set(gca,  'units', 'pixels', 'innerposition', [169.666666666667,100,736.333333333333,200])
end



% Histograms
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(periodD * 60, 'binwidth', 8, 'normalization', 'probability');
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0.0; h.LineWidth = 1.5;
xlabel('\tau_{\itd} (min)')
xt = 0:120:360; xticks(xt)
xlim([0 360])
scatter(nanmean(periodD(:)*60), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 470, 4, 40)

new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(pcorrel, 'binwidth', 0.02, 'normalization', 'probability');
h.FaceColor = [0.00,0.45,0.74]; h.EdgeAlpha = 0.0; h.LineWidth = 1.5;
xlabel('\rho')
xlim([-0.45 0.45])
scatter(nanmean(pcorrel), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 440, 4, 40)

%% Quad order



%% Scatter (save)
figure;
L = 20; i = 1;
for t = t0 : (t0+dur)
% for t = t0
    x = qs.pos(:, cell_ind, 1, t);
    y = qs.pos(:, cell_ind, 2, t);
    mx = nanmean(x); my = nanmean(y);
    scatter(x, y, 50, 'markerfacecolor', 'k')
    for j = 1:4
        text(x(j), y(j), num2str(j), 'fontsize', 20, 'color', 'r')
    end
    xlim([mx-L mx+L])
    ylim([my-L my+L])
    saveas(gca, ['img_' num2str(i) '.png'])
    i = i + 1;
end

%% Load and adjust CPM image sequence
% Read the image stack
IMG = [];
for i = 1:30
    img = imread("C:\Users\HyunGyu\mda\CellCluster Class\simulation trajectories\S 2.8 E -65\quad_cpm(3400_to_3700).tif", i);
    img = reshape(img, [1 size(img)]);
    IMG = cat(1, IMG, img);
end

% figure; imshow(squeeze(IMG(109, :, :, :)));
white = [255; 255; 255];
black = [0; 0; 0];
yellow = [255; 255; 0]; % 1
blue = [0; 0; 255]; % 3
green = [0; 100; 0]; % 2
orange = [255; 165; 0]; % 4

% Remove the black line & fill it with color
IMG(:, 22, :, :) = 255;
for t = 1:size(IMG,1)
    for i = 2:size(IMG,2)-1
        a = squeeze(IMG(t, 22, i, :));
        fin = 0; flag0 = 0; flag1 = 0;
        for x = [-1 0 1]
            for y = [-1 0 1]
                if ~(x==0 && y==0) && ~(x==1 && y==0)
                    an = squeeze(IMG(t, 22+y, i+x, :));
                    if isequal(an, white)
                        flag0 = 1;
                    end
                    if isequal(an, blue) || isequal(an, orange) || isequal(an, green) || isequal(an, yellow)
                        flag1 = 1;
                    end
                    if flag0 == 1 && flag1 == 1
                        fin = 1;
                        break;
                    end
                end
            end
            if fin == 1
                IMG(t, 22, i, :) = black;
                break;
            end
        end
        if flag0 == 0 && flag1 == 1
            IMG(t, 22, i, :) = IMG(t, 23, i, :);
        end
    end
end

%% Using 'remove_boundary' function
% Fresh import IMG from .tif
wb = waitbar(0);
for t = 1:30
    img = squeeze(IMG(t, :, :, :));
    IMG(t, :, :, :) = remove_boundary(qs, img, 1);
    waitbar(t/30)
end
close(wb)

%%
% Change color to match (blu, ora, yel, pur)
factor = 1.5;
% factor = 0.6;
for t = 1:size(IMG,1)
for i = 1:size(IMG, 2)
    for j = 1:size(IMG, 3)
        a = squeeze(IMG(t, i, j, :));
        if ~isequal(a, white) && ~isequal(a, black)
            if isequal(a, yellow) 
                IMG(t, i, j, :) = min(qs.blu * 255 * factor, 255);
            elseif isequal(a, green) 
                IMG(t, i, j, :) = min(qs.ora * 255 * factor, 255);
            elseif isequal(a, blue) 
                IMG(t, i, j, :) = min(qs.yel * 255 * factor, 255);
            elseif isequal(a, orange) 
                IMG(t, i, j, :) = min(qs.pur * 255 * factor, 255);
            else
                IMG(t, i, j, :) = [50 50 50];
            end
        elseif isequal(a, white)
            IMG(t, i, j, :) = black;
        end
    end
end
end
              
% Save the image sequence
for t = 1:30
%     imshow(squeeze(IMG(t, :, :, :)));
%     saveas(gca,['img_' num2str(t) '.png'])
    imwrite(squeeze(IMG(t, :, :, :)), ['img_' num2str(t) '.png']); % Save matrix as an image directly
end

%% One snapshot (beside 3D tube)
img = squeeze(IMG(1, :, :, :));
factor = 1.5;

blu = uint8(min(qs.blu * 255 * factor, 255));
ora = uint8(min(qs.ora * 255 * factor, 255));
yel = uint8(min(qs.yel * 255 * factor, 255));
pur = uint8(min(qs.pur * 255 * factor, 255));

for i = 2:size(img, 1)-1
    for j = 2:size(img, 2)-1
        a = squeeze(img(i, j, :));
        if ~isequal(a, white) && ~isequal(a, black)
            if isequal(a, yellow) 
                img(i, j, :) = min(qs.blu * 255 * factor, 255);
            elseif isequal(a, green) 
                img(i, j, :) = min(qs.ora * 255 * factor, 255);
            elseif isequal(a, blue) 
                img(i, j, :) = min(qs.yel * 255 * factor, 255);
            elseif isequal(a, orange) 
                img(i, j, :) = min(qs.pur * 255 * factor, 255);
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
            img(i,j,:) = n1;
        end
    end
end

% Remove stupid colors 
img_copy = img;
for i = 2:size(img, 1)-1
    for j = 2:size(img, 2)-1
        a = squeeze(img_copy(i, j, :)); a = a';
        if ~isequal(a, white') && ~isequal(a, blu) && ~isequal(a, ora) && ~isequal(a, yel) && ~isequal(a,pur)
            b = img(i-1:i+1, j-1:j+1, :);
            b = reshape(b, 9, 3);
            b(5,:) = [];
            
            wt = sum(b == white', 2);
            bl = sum(b == blu, 2);
            or = sum(b == ora, 2);
            ye = sum(b == yel, 2);
            pu = sum(b == pur, 2);
            
            wt_n = sum(wt == 3);
            bl_n = sum(bl == 3);
            or_n = sum(or == 3);
            ye_n = sum(ye == 3);
            pu_n = sum(pu == 3);
            
            [~, ind] = max([wt_n, bl_n, or_n, ye_n, pu_n]);
            
            if ind == 1
                img(i,j,:) = white;
            elseif ind == 2
                img(i,j,:) = blu;
            elseif ind == 3
                img(i,j,:) = ora;
            elseif ind == 4
                img(i,j,:) = yel;
            elseif ind == 5
                img(i,j,:) = pur;
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

%% Smooth curves
% scale = 250 / 120;
scale = 2.6;
quad_s = qs.pos(1:4, cell_ind, :, t0:(t0+dur));
quad_s = permute(quad_s, [4,3,2,1]);
quad_s = squeeze(quad_s) * scale;
offset = nanmean(quad_s, 3);
quad_s = quad_s - offset(1, :); 
quad_s(:,2,:) = - quad_s(:,2,:); % Invert y
% quad_s(:,1,:) = quad_s(:,1,:) + shiftx;
% quad_s(:,2,:) = quad_s(:,2,:) + shifty;

[T, ~, N] = size(quad_s);
q_smooth_s = nan(size(quad_s));
for i = 1:4
    for j = 1:2
        y = quad_s(:, j, i); ind = find(~isnan(y));
        x = 1:length(y);
        y = y(ind); x = x(ind);
%         f = fit(x', y, 'smoothingspline', 'smoothingparam', 0.01);
        f = fit(x', y, 'smoothingspline', 'smoothingparam', 0.06);
        q_smooth_s(ind, j, i) = f(x);
    end
end

    
%% Overlay trajectory
ms = 300; lw = 10; malpha = 0.8;
L = 35;

% old 
% t00 = 121; % 4 hour (mix)
% t00 = 210; % 7 hour (rigid)

% t00 = 1; % nomix (0 hour)
t00 = 110; % mix (3.667 hour)
for dt = [0 30 60 90]
    figure;
    t = t00 + dt;
    img = squeeze(IMG(floor(t/10)+1, :, :, :));
    % l1 = round(qcm(t,:) - L); l2 = round(qcm(t,:) + L);
    center = round(cpos(floor(t/10)+1, :));
    shiftx = 87 - center(1);
    shifty = 87 - center(2);
    img = img(center(2)-L:center(2)+L, center(1)-L:center(1)+L, :);
    imshow(img);
    hold on
    for i = 1:4
        x = q_smooth_s(t00:t, 1, i) + shiftx + 2;
        y = q_smooth_s(t00:t, 2, i) + shifty + 1;
        plot(x, y, 'linewidth', lw, 'color', [qs.colors(i, :) malpha]);
        scatter(x(end), y(end), ms, 'markerfacecolor', qs.colors(i, :), 'markeredgecolor', 'none');
    end
    set(gca, 'position', [0,0,1.0,1.0])
    set(gcf, 'pos', [500 500 300 300])
    % saveas(gcf, ['fig_' num2str(dt) '(mix).png'])
    % saveas(gcf, ['fig_' num2str(dt) '(nomix).png'])
end

%% Get center from image
x = nan(4, 100); y = nan(4, 100);
n = ones(4, 1);
for i = 1:100
    for j = 1:100
        a = squeeze(IMG(t, i, j, :)); a = a';
        if isequal(a, uint8(qs.blu * 255))
            x(1, n(1)) = j;
            y(1, n(1)) = i;
            n(1) = n(1) + 1;
        elseif isequal(a, uint8(qs.ora * 255))
            x(2, n(2)) = j;
            y(2, n(2)) = i;
            n(2) = n(2) + 1;
        elseif isequal(a, uint8(qs.yel * 255))
            x(3, n(3)) = j;
            y(3, n(3)) = i;
            n(3) = n(3) + 1;
        elseif isequal(a, uint8(qs.pur * 255))
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

%% Get center for the all cells (non-white (or black) regions)
cpos = nan(30, 2);
for time = 1:30
x = []; y = [];
for i = 1:100
    for j = 1:100
        a = squeeze(IMG(time, i, j, :));
%         if ~isequal(a, uint8(white))
        if ~isequal(a, uint8(black))
            x = [x j];
            y = [y i];
        end
    end
end
cpos(time,1) = mean(x);
cpos(time,2) = mean(y);
waitbar(time/301)
end

%% Plot 3d
cell_ind = 2; 
% t0 = 1800; dur = 300; 
t0 = 1; dur = 3000;
timespan = t0 : (t0+dur);
rect = [1400 200 700 850];
lim = 50;
viewangle = [327.6,379.0,41.9];
figure('pos', rect); hold on 
plot3d_traj(qs, cell_ind, timespan, rect, lim, viewangle, []);
% xticks(-50:25:50);
% yticks(-50:25:50);
% xticklabels(["" "" ""])
% yticklabels(["" "" ""])
xticks(-15:15:15);
yticks(-15:15:15);
img = squeeze(IMG(1, :, :, :));
% img = permute(img, [2 1 3]);
img = img(end:-1:1, :, :);
[X,Y] = meshgrid(-49:50, -49:50); Z = zeros(100, 100); 
X = X / (250 / 120); Y = Y / (250 / 120);
X = X - 1;
Y = Y - 1;
s = surf(X, Y, Z, img); s.EdgeColor = 'none';
% xlabel('x')
% ylabel('y')

%% 3d Tube
cell_ind = 2; 
t0 = 2400; dur = 300;
timespan = t0 : (t0+dur);
rect = [0 0 700 850];
lim = 16;
viewangle = [221.2025112125169,246.8015159319472,90.59189229950624];
new_fig(rect, 4.5, 45, [], [], [], [], "");
interpN = 4; 
twfactor = 7; 
zfactor = 0.1;
sp = 0.008;
tube_surf(qs, cell_ind, timespan, rect, lim, viewangle, twfactor, interpN, sp, zfactor);
% Overlay bottom
img_ground = zeros(200,200,3); Z0 = zeros(200,200); 
[X,Y] = meshgrid(-99:100, -99:100);
s = surf(X,Y,Z0, img_ground); s.EdgeColor = 'none';
img = squeeze(IMG(1, :, :, :));
img = img(end:-1:1, :, :);
[X,Y] = meshgrid(-49:50, -49:50); Z = zeros(100, 100) + 0.01; 
X = X / scale; Y = Y / scale;
X = X - 1;
Y = Y + 0;
s = surf(X, Y, Z, img); s.EdgeColor = 'none';
xt = -40/3:40/6:40/3; xticks(xt)
yt = -40/3:40/6:40/3; yticks(yt)
zticks(0:7.5:30);
zticklabels(string(0:150:600))
xticklabels(repmat("", [3 1]))
yticklabels(repmat("", [3 1]))
zlim([0 10*30*zfactor])

f = fill3([-lim lim lim -lim], [-lim -lim -lim -lim]+0.01, [0 0 3 3]*60/20, ...
    [0.07,0.62,1.00]); f.FaceAlpha = 0.2; f.EdgeColor = 'none';
f = fill3([-lim -lim -lim -lim]+0.01, [-lim lim lim -lim], [0 0 3 3]*60/20, ...
    [0.07,0.62,1.00]); f.FaceAlpha = 0.2; f.EdgeColor = 'none';

f = fill3([-lim lim lim -lim], [-lim -lim -lim -lim]+0.01, [3.66 3.66 6.66 6.66]*60/20, ...
    [1.00,0.41,0.16]); f.FaceAlpha = 0.2; f.EdgeColor = 'none';
f = fill3([-lim -lim -lim -lim]+0.01, [-lim lim lim -lim], [3.66 3.66 6.66 6.66]*60/20, ...
    [1.00,0.41,0.16]); f.FaceAlpha = 0.2; f.EdgeColor = 'none';

camroll(-90)
setfigprop(1300,1200,4,40)
ztickangle(-90)
set(gca, 'innerposition', [0.13,0.11,0.775,0.74], 'position', [0.13,0.11,0.775,0.74])
zl = zlabel('');



%% Angle (all) and period histo
a = squeeze(nanmean(th, 1));
figure('pos', [293,609,694,500]); hold on 
% figure('pos', [100,100,694,600]); hold on 
lw = 1.5; alpha = 0.5;
for i = 1:100
    x = (0:450-1) / 30;
    y = a(i, 1:450);
    col = rand(1,3);
    plot(x, y, 'linewidth', lw, 'color', [col alpha]); 
end
set(gca, 'linewidth', 3.7, 'fontsize', 37)
yt = (-10:5:10)*pi; yticks(yt);
ytl = string(-10:5:10) + "\pi"; ytl(3) = "0"; yticklabels(ytl);
ylim([-10*pi 10*pi])
xlim([0 15])
xlabel('time')
ylabel('\theta (rad)')

% Get 'A'
Dt = 1:10:401;
A = nan(length(Dt), qs.N, qs.T);
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

%% Quadorder & transition
% Import 'quadorder_simul.csv' as table
quadordersimul = quadordersimul{:, :};
quadordersimul = quadordersimul(2:end,:);
figure;
for i = 1:qs.N
    subplot(10, 10, i)
    plot(quadordersimul(i, :))
    title(i)
end

new_fig([50 50 500 500], 1.5, 15, [], [], [], [], "");
h = histogram(quadordersimul(:), 'normalization', 'probability');
xticks(1:6); xlabel("state")
dt = 10;
trans = [1 2; 1 3; 1 4; 1 5; 1 6; 2 3; 2 4; 2 5; 2 6; 3 4; 3 5; 3 6; 4 5; 4 6; 5 6];
quadorder_trans = nan(qs.N, size(quadordersimul, 2)-dt);
for n = 1:qs.N
    y = squeeze(quadordersimul(n, :));
    for t = 1:size(quadordersimul, 2)-dt
        qo0 = y(t);
        qo1 = y(t+dt);
        if qo0 ~= qo1
            for i = 1:6
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

qos = quadorder_trans(:);
ind = 1:15; ind1 = [1 10 15];
ind = setdiff(ind, ind1);
for i = 1:length(qos(:))
    if ismember(qos(i), ind)
        qos(i) = 1;
    elseif ismember(qos(i), ind1)
        qos(i) = 2;
    end
end
qos(qos == 0) = [];
qos(isnan(qos)) = [];

new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
% h = histogram(qos, ...
%     'normalization', 'probability');
% h.FaceColor = [1.00,0.41,0.16]; h.EdgeAlpha = 0.1; h.LineWidth = 2.5;
%  xlabel('transition type')
h = histcounts(qos); h = h / sum(h);
b = bar(1:2, h);
b.BarWidth = 0.4; b.FaceAlpha = 0.5; b.EdgeColor = 'none';
xticks([1,2]);
xlim([0.5 2.5])
xticklabels({"", ""})
setfigprop(500, 328, 4, 40)
        
trans_freq = sum(quadorder_trans ~= 0 & ~isnan(quadorder_trans), 2) ... 
    ./ sum(~isnan(quadorder_trans), 2) * 30;

new_fig([50 50 600 600], 1.5, 15, [], [], [], [], "");
histogram(trans_freq, 'binwidth', 0.35, ...
    'normalization', 'probability'); xlabel("Mean  # of exchange per hour")


%% Unwrapping the theta vs. time
th = squeeze(nanmean(th, 1));
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
 t_turn_q = nan(100,200); 
 t_turn_q(1, :) = 1;
 for i = 1:200
     y = th(i, :);
     
     N = 1;
     t0 = t_turn_q(1, i);
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
         t_turn_q(N, i) = t_new;
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
tt = t_turn_q(:, i); tt(isnan(tt)) = [];
plot(time, th(i, :), 'linewidth', 2); 
hold on
plot(time(tt), th(i, tt), 'x', 'linewidth', 2)

%% Get the slopes
w = nan(100, 200);
for i = 1:200
    a = t_turn_q(:, i);
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
% mw = nanmean(abs(w), 1);
% figure; histogram(mw)
% xlabel('<|\omega|>')
% 
% mperiod = 2*pi ./ mw;
% figure; histogram(mperiod*2)
% xlabel('2\pi/|\omega|')        
% xlim([0 720])

mw = nanmean(abs(w), 1);
mperiod = 2*pi ./ mw;

new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(mperiod*2, 'binwidth', 8, 'normalization', 'probability');
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0.0; h.LineWidth = 1.5;
xlabel('{\itT} (min)')
xticks(0:240:720)
% xlim([50 1500])
xlim([0 720])
scatter(nanmean(mperiod*2), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 470, 4, 40)

%% Unwrapping the whole theta (based on initial slope)
th_unwrap = th;
for i = 1:200
    sig = sign(w(1, i));
    a = t_turn_q(:, i);
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
ylabel('<\theta> (rad)')
% ylim([-16*pi 16*pi])
ylim([-12*pi 12*pi])
yt = [-12*pi : 6*pi : 12*pi]; yticks(yt);
% ytl = {'-16\pi' '-8\pi' '0' '8\pi' '16\pi'}; yticklabels(ytl)
yt = string(-12:6:12) + "\pi"; yt(3) = '0'; yticklabels(yt)
xt = 0:300:1200; xticks(xt)
box off
setfigprop(1200, 660, 4, 40)

% figure('pos', [0,0,693,479]); hold on 
% for c=1:40
%     plot((0:558-1)/30, thFlat(c, :), 'linewidth', 2.4, 'color', [rand(1,3) 0.6]);
% end
% xlim([0 20])
% r = refline(0, 0); r.Color = 'k'; r.LineWidth = 3; r.LineStyle = '--';
% xlabel('t (min)')
% ylabel('\theta (rad)')
% ylim([-60 60])
% set(gca, 'fontsize', 37, 'linewidth', 3.7)
% yt = [-16*pi : 4*pi : 16*pi]; yticks(yt);
% ytl = {'-16\pi' '-12\pi' '-8\pi' '-4\pi' '0' '4\pi' '8\pi' '12\pi' '16\pi'}; yticklabels(ytl)
% xt = 0:10:20; xticks(xt)
% xt = xt * 60; xticklabels(xt)
% ylim([-8*pi 8*pi])
% 
% box off




