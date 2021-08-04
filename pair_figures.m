% Pair figures (experiments)
p = CellCluster(pos_pair);

dth = ang_vel(p);
D = cell_dist(p);
th = ang(p);
dt = 5;
dth_2 = nan(2, 57, p.T-dt);
for t = 1:p.T-dt
    dth_2(:, :, t) = th(:, :, t+dt) - th(:, :, t);
end
dth_2 = dth_2 / dt;

pc = pair_correlation(p);
mpc = nanmean(pc, 2);

newD = subtract_smoothing(p, D, 3e-6);
newDth = subtract_smoothing(p, dth, 3e-6) * 30; % rad per hour

periodD = get_period_fft(p, newD); % In hour
periodDth = get_period_fft(p, newDth); % hour

pcorrel = pearson_correlation(p, newD, abs(newDth));

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
plot3d_traj(p, cell_ind, timespan, rect, lim, viewangle, []);
img = imread('..\COMBINED_pair_Oct8_Oct23_Oct15(59).tif', 1);
hold on; overlay_image(p, img, cell_ind, lim, 0, 4);

% 3D Tube
cell_ind = 31; 
timespan = 1:668;
rect = [0 0 700 850];
lim = 30;
viewangle = [389.1,330.9,211.3];
new_fig(rect, 4.5, 45, [], [], [], [], "");
interpN = 4; twfactor = 5; sp = 0.008;
zfactor = 0.2;
tube_surf(p, cell_ind, timespan, rect, lim, viewangle, twfactor, interpN, sp, zfactor);
img = imread('..\COMBINED_pair_Oct8_Oct23_Oct15(59).tif', 1);
hold on; overlay_image(p, img, cell_ind, lim, 0, 4);
xticklabels(repmat("", [3 1]))
yticklabels(repmat("", [3 1]))
zlim([0 15*30*zfactor])
camroll(-90)
setfigprop(1700,1200,4,40)
ztickangle(-90)
% Mark reverse-turn
% scatter3(15, -lim, 24, 1000, '>', 'r', 'filled')
% tubeplot(x2, y2, z, 0.5, ones(length(z),1), 10);


%%
% dth_smooth = smoothing(p, dth, 0.45);
dth_2_smooth = smoothing(p, dth_2, 0.85);
D_smooth = smoothing(p, D, 0.1);
th_smooth = smoothing(p, th, 0.75);

% theta vs. time
time = (0:p.T-1)/30;
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], "");
% yyaxis left; 
ylabel('\theta (rad)')
plot(time, squeeze(th_smooth(1, cell_ind, :) + pi/6), 'linewidth', 6, 'color', p.colors(1,:));
ylim([-6*pi 6*pi])
yt = -3*pi : 3*pi : 3*pi; yticks(yt);
ytl = {'-3\pi', '0', '3\pi'}; yticklabels(ytl);
% yyaxis right; ylabel('\omega (rad/h)'); ylim([-10 10])
% plot(time(1:end-dt), squeeze(dth_2_smooth(1, cell_ind, :)*30), 'linewidth', 3, 'color', [0.3 0.3 0.3 0.9]);
% plot(time(1:end-dt), squeeze(dth_2(1, cell_ind, :))*30, '.', 'markersize', 15, 'color', [0.3 0.3 0.3 0.5]);
r = refline(0,0); r.LineWidth = 3; r.LineStyle = '--'; r.Color = p.gra;
yt = -4*pi : 4*pi : 5*pi; yticks(yt);
yt = string(-4:4:4) + '\pi'; yt(2) = '0'; yticklabels(yt);
xt = 0:5:15; xticks(xt);
xt = xt * 60; xticklabels(xt)
% xlabel('t (min)')
box on
% Mark reverse-turn
scatter(4, 10, 500, '^', 'k', 'filled')
setfigprop(1000, 450, 4, 40)

% omega vs. D
ms = 100; 
lw = 4;
% lw = 6;
time = (0:p.T-1)/30;
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], "");
yyaxis left; ylabel('|\omega| (rad/min)')
plot(time(1:end-dt), squeeze(abs(dth_2_smooth(1, cell_ind, :)) * 30 / 60), 'linewidth', lw,  'color', [p.gra 0.9]);
scatter(time(1:end-dt), squeeze(abs(dth_2(1, cell_ind, :)) * 30 / 60), ms,  'markerfacecolor', p.gra, ...
    'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', p.gra)
% yt = 0 : 2*pi : 4*pi; yticks(yt);
% ytl = string(0:2:4) + "\pi"; ytl(1) = '0'; yticklabels(ytl)
ylim([0 0.23])
% ylim([0 2*pi+1.0])
yyaxis right; ylabel('\itd (\mum)')
plot(time, squeeze(D_smooth(1, cell_ind, :)), 'linewidth', lw, 'color', p.gre);
scatter(time, squeeze(D(1, cell_ind, :)), ms , 'markerfacecolor', p.gre, ...
        'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', p.gre)
% xlabel('t (hour)')
yt = 6:4:14; yticks(yt);
xt = 0:5:15; xticks(xt); xt = string(xt * 60); xticklabels(xt) 
% xlim([5 10]); 
% ylim([6 13])
box on
setfigprop(1040, 450, 4, 40)
set(gca, 'innerpos', [0.15,0.2502,0.711858981952835,0.749814818957713])

% omega vs. D (inset)
ms = 200; 
lw = 8;
time = (0:p.T-1)/30;
new_fig([50 50 1000 500], 4.5, 45, [], [], [], [], "");
yyaxis left; ylabel('|\omega| (rad/min)')
plot(time(1:end-dt), squeeze(abs(dth_2_smooth(1, cell_ind, :)) * 30 / 60), 'linewidth', lw,  'color', [p.gra 0.9]);
scatter(time(1:end-dt), squeeze(abs(dth_2(1, cell_ind, :)) * 30 / 60), ms,  'markerfacecolor', p.gra, ...
    'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', p.gra)
% yt = 0 : pi : 4*pi; yticks(yt);
% ytl = {'0', '\pi', '2\pi', '3\pi', '4\pi'}; yticklabels(ytl);
% ylim([0 4*pi])
% ylim([0 0.23])
% ylim([0 2*pi+1.0])
yyaxis right; ylabel('\itd (\mum)')
plot(time, squeeze(D_smooth(1, cell_ind, :)), 'linewidth', lw, 'color', p.gre);
scatter(time, squeeze(D(1, cell_ind, :)), ms , 'markerfacecolor', p.gre, ...
        'markeredgecolor', 'none', 'markerfacealpha' , 0.5);
set(gca, 'YColor', p.gre)
% xlabel('t (hour)')
% yt = 6:4:14; yticks(yt);
yt = 6:3:15; yticks(yt);
% xt = 0:5:15; xticks(xt); xt = string(xt * 60); xticklabels(xt) 
xt = 5:2.5:10; xticks(xt); xt = string(xt * 60); xticklabels(xt)
xlim([5 10])
% ylim([6 13])
setfigprop(1040, 400, 4, 40)

% Histograms
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(periodD*60, 'binwidth', 30, 'normalization', 'probability');
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0.1; h.LineWidth = 1.5;
xlabel('\tau_{\itd} (min)')
xt = 0:120:480; xticks(xt);
% xlim([0 390])
scatter(nanmean(periodD(:)*60), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 470, 4, 40)

new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(pcorrel, 'binwidth', 0.07, 'normalization', 'probability');
h.FaceColor = [0.30,0.75,0.93]; h.EdgeAlpha = 0.1; h.LineWidth = 1.5;
xlabel('\rho')
xlim([-0.6 0.6])
scatter(nanmean(pcorrel(:)), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 440, 4, 40)

% Histogram of pair correlation
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(mpc, 'binwidth', 0.08, 'normalization', 'probability');
h.FaceColor = [1.00,0.41,0.16]; h.EdgeAlpha = 0.1; h.LineWidth = 1.5;
xlabel('\it{c}_{pair}')
xlim([-0.7 0.7])
scatter(nanmean(mpc(:)), 0.0, 350, [1 0 0], '^', 'filled')
setfigprop(500, 470, 4, 40)


%%
figHandles = findall(0,'Type','figure'); 

savefig(figHandles, "figures\figs_pair_dec_14.fig");



























%% 0. Calculate the angle
position = pos_pair;
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
dth = squeeze(dth(1, :, :));
dth = dth;

for t=1:T-1
    u1 = du(:,:,:,t);
    u2 = du(:,:,:,t+1);
    crossu = u1(:,:,1).*u2(:,:,2) - u1(:,:,2).*u2(:,:,1);
    crossu = - crossu;
    ang = asin(crossu);
    dth(:,t) = ang(1,:);
%     waitbar(t/T)
end

sc = 1000 / 340.5;
D = position(1, :, :, :) - position(2, :, :, :);
D = absvec(D, 3) / 2;
D = squeeze(D) * sc; % Calibrated

sz = size(dth);
sz(2) = sz(2)+1;
th = nan(sz);
th(:,1) = 0;

for i=1:size(position,2)
    a = dth(i,:);
    a = squeeze(a);
    a = a(~isnan(a));
    ac = cumsum(a);
    th(i,1:length(ac)) = ac; % Start time doesn't matter much
end

%% 1. FFT
%% Subtract heavily smoothed D
newD = D;
for i=1:size(D, 1)
    y = D(i,:); y = y';
    IND = find(~isnan(y));
    x = 1:length(y); x = x';
    x = x(IND);
    y = y(IND);
    if ~isempty(y)
        f = fit(x, y, 'smoothingspline', 'SmoothingParam', 3e-6);
        newD(i, IND) = newD(i, IND) - reshape(f(x), 1, length(IND));
    end
end

% Check 
figure; 
for i=1:size(newD, 1)
    subplot(10, 6, i)
    plot(newD(i,:))
end
%% Subtract heavily smoothed dth
newDth = dth;
for i=1:size(dth, 1)
    y = dth(i,:); y = y';
    IND = find(~isnan(y));
    x = 1:length(y); x = x';
    x = x(IND);
    y = y(IND);
    if ~isempty(y)
        f = fit(x, y, 'smoothingspline', 'SmoothingParam', 3e-6);
        newDth(i, IND) = newDth(i, IND) - reshape(f(x), 1, length(IND));
    end
end

% Check 
figure; 
for i=1:size(newDth, 1)
    subplot(10, 6, i)
    plot(newDth(i,:))
end

%% FFT subplots
figure;
N = size(D, 1);
T = 2/60; % Interval (hours)
Fs = 1/T;
maxT = nan(N ,1);
j = 1;
for i=1:1:57
    
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

        subplot(10, 6, j)
        
        yyaxis left
        plot(t, P1, '.-', 'markersize', 10, 'linewidth', 2)
%          
        % Smooth the curve
        f = fit(t(2:end)', P1(2:end)', 'smoothingspline', 'SmoothingParam', 0.90);
        f = f(t(2:end));
        ind = find(max(f) == f) + 1;
        hold on 
        plot(t(2:end), f, 'g')
        plot(t(ind), 0, 'kx', 'linewidth', 2, 'markersize', 8)
        maxT(i) = t(ind);
        title(['T_r_g= ' num2str(round(t(ind),3))])



        yyaxis right
        plot((0:L-1)*T, dd)
        ylabel('r_g')


   
        xlim([0 15])
        
        j = j + 1;
    end
    
end


%% Display the period 
figure('pos', rect); 
h = histogram(maxT, 'binwidth', 0.5);
h.FaceColor = 'k'; h.FaceAlpha = 0.7; h.EdgeAlpha = 0;
% title('Oscillation period of r_g')
xlabel('T (hour)')
ylabel('n')
set(gca, 'fontsize', 35, 'linewidth', 3.5)
box off


%% 2. The angles vs time
T = size(th, 2);
figure; hold on
ii = 1;
% q = [11, 30, 25];
% p = [29 43 33];
t = 1:668;
for pind = 31
% for i=1:10
    fill([5 10 10 5]*30, [-30 -30 30 30], 'r')
%     subplot(3, 1, ii)
%     subplot(7,9,ii)
%         plot((0:T-1)*2/60, squeeze(th(j,pi,:)), 'linewidth', 4, 'color', 'k'); hold on 
%     plot((0:T-1)*2/60, squeeze(th(pi,:)), 'linewidth', 4); hold on 
    
    y1 = -th(pind,t) + th(pind, t(1)) + pi/2;
    y2 = -th(pind,t) + th(pind, t(1)) + pi/2 - pi;

    
    f1 = fit(t(~isnan(y1))', y1(~isnan(y1))', 'smoothingspline', 'smoothingparam', 0.005);
    f2 = fit(t(~isnan(y2))', y2(~isnan(y2))', 'smoothingspline', 'smoothingparam', 0.005);
    
    t = t(~isnan(y1));
    y1 = f1(t);
    y2 = f2(t);
    
    plot(t, y1, 'linewidth', 4); hold on 
    plot(t, y2, 'linewidth', 4);
    
%     plot((1:T), -th(pind,:), 'linewidth', 2); hold on 
%     plot(1:T, -(th(pind, :) - pi), 'linewidth', 2);
%     ylim([-40 40])
    
%         ylim([-15 15])
%     xlim([2 14])
%     xticks = 2:2:14;
%     xticklabels(string(0:2:12))
%     r = refline(0, 0); r.LineStyle = ':'; r.Color = [0 0 0];
%     set(gca, 'fontsize', 25, 'linewidth', 3)
%     box off
%     ylabel('\theta (rad)')
    ii = ii + 1;
%     title(num2str(pind))
%     refline(0, 0)
end
xt = t(1):120:t(end);
xtl = string((0:120:800)/30);
xticks(xt);
xticklabels(xtl);


ylim([-6*pi 6*pi])

yt = [-6*pi : 2*pi : 6*pi]; yticks(yt);
ytl = {'-6\pi' '-4\pi' '-2\pi' '0' '2\pi' '4\pi' '6\pi'}; yticklabels(ytl)

set(gca, 'fontsize', 35, 'linewidth', 3.5) 
ylabel('\theta')
xlabel('hour')

%% Step 1.
% Manually find locations where slopes change significantly
M = {};
M{1} = [1 125 184];
M{2} = [1 33 96 127 277 298 413 447 504 544 595 656];
M{3} = [1 113 146 158 291 460];
M{4} = [1 32 92 169 203];
M{5} = [15 121 204 232 391];
M{6} = [1 24 42 123 150 330 342 398 409 447];
M{7} = [1 37 65 125 200 295 317];
M{8} = [1 35 103 140 210 240];
M{9} = [1 29 105 113 170];
M{10} = [1 40 58 139 168 375 413 447 470 515 604];
M{11} = [1 136 203 234 365 429];
M{12} = [1 43 145 176 191 253 269];
M{13} = [1 59 74 88 99 121 125 140 158 283];
M{14} = [1 111 180 226 237 297 316 336 373 297];
M{16} = [1 16 52 64 90 110 123 158 195 255 267 289 306 315 432]; 
M{17} = [1 72 124 227 392];
M{18} = [1 73 104 113 250 373];
M{19} = [1 27 82 178];
M{20} = [1 142 341];
M{21} = [1 33 59 140 173 264 306 353 373];
M{22} = [1 54 125 382];
M{23} = [1 104 216 248 396];
M{24} = [1 127 313 443];
M{25} = [1 267 456];
M{26} = [1 81 170 269 390];
M{27} = [1 48 63 136 164 193 228 246 278 341 556];
M{28} = [1 55 124 165];
M{29} = [1 33 86 339 483];
M{30} = [1 30 103 280 346];
M{31} = [1 124 442];
M{32} = [1 161 322 455];
M{33} = [1 223 529];
M{34} = [1 21 86 138 152 190 230 274 298 383 447];
M{35} = [1 63 95 201 228];
M{36} = [1 16 49 62 67 105];
M{37} = [1 27 90 191];
M{38} = [1 55 110 185 197 271 283 335 354 388];
M{39} = [1 55 81 158 198 345];
M{40} = [1 60 121 148 234 258 311];
M{41} = [1 422];
M{42} = [1 45 166 243 385];
M{43} = [1 55 71 175 191 246 255 556];
M{44} = [1 23 84 164 272];
M{45} = [1 16 36 53 92 106 120 130 136 191 203 226 244 253 299 309 329 357 389 394 464 486];
M{46} = [1 413];
M{47} = [1 121 161 174 201 216 252  454 533];
M{48} = [1 49 63 84 101 144 302 344 392];
M{49} = [1 120 196 334 433];
M{50} = [1 204];
M{51} = [1 30 54 62 88 118 129 226];
M{52} = [1 352 437];
M{53} = [1 31 50 76 122 185 243];
M{54} = [1 52 80 108 127 147 154 171 219 281 324 345 423];
M{55} = [1 13 27 34 47 66 118 131 188 215 224 258 266 281 311 324 355 433];
M{56} = [1 13 22 50 60 98 105 148 159];
M{57} = [1 54 80 94 110 134 167 198];

th = squeeze(th(1, :, :));
dth = squeeze(dth(1, :, :));
% Start ~ end if not specified
for c = 1:57
    if isempty(M{c})
        i = sum(~isnan(th(c,:)));
        M{c} = [1 i];
    end
end

% Visualize
T = size(th, 2);
figure; hold on
ii = 1;
for pind = 1:57
    subplot(7,9,ii)
    plot((0:T-1)*2, squeeze(th(pind,:)), 'linewidth', 2); hold on 
    if ~isempty(M{pind})
        for m = M{pind}
            plot(m * 2, th(pind, m), 'x', 'linewidth', 2, 'markersize', 6)
        end
    end
    ylim([-40 40])
    ii = ii + 1;
%     title(num2str(pind))
    r = refline(0, 0); r.LineWidth = 1; r.LineStyle = '--';
    xlabel('min')
    ylabel('\theta')
end

%% Step 2 : Find slopes between each specified points
slope = {};
for c = 1:57
    dy = diff(th(c, M{c})); % In radians
    dx = diff(M{c})/30; % in hours
    slope{c} = dy./dx;
end

% Concatenate all the slopes
allSlopes = [];
for c = 1:57
    allSlopes = cat(2, allSlopes, slope{c});
end

% Display the histogram
figure; histogram(abs(allSlopes(:)))
figure; 
for i=1:320
    refline(allSlopes(i), 0)
end

%% Step 3 : For each cell, remove similar consecutive points (redunant) 
% whose difference of slope is smaller than theshold

Mnew = M;
threshDiff = 0.7;
for c=1:57
    m = M{c};
    diffm = diff(th(c, m), 2);
    idx = [];
    for i = 2:length(m)-1
        if abs(diffm(i-1)) <= threshDiff
            idx = [idx i];
        end
    end
    m(idx) = [];
    Mnew{c} = m;
end

% Display the result
figure; hold on
ii = 1;
for pind = 1:57
    subplot(7,9,ii)
    plot((1:T), squeeze(th(pind,:)), 'linewidth', 2); hold on 
    if ~isempty(Mnew{pind})
        % Original pts
        for m = M{pind}
            plot(m, th(pind, m), 'bx', 'linewidth', 2, 'markersize', 6)
        end
        
        % Modified pts
        for m = Mnew{pind}
            plot(m, th(pind, m), 'r.', 'linewidth', 2, 'markersize', 15)
        end
    
    end
    ylim([-40 40])
    ii = ii + 1;
    title(num2str(pind))
    refline(0, 0)
end

%% Step 4. Slopes of Mnew
slope = {};
for c = 1:57
    dy = diff(th(c, Mnew{c})); % In radians
    dx = diff(Mnew{c})/30; % in hours
    slope{c} = dy./dx;
end

slope_change_freq = [];
for i = 1:57
    a = slope{i};
    a = sign(a);
    b = diff(a);
    b = sum(b ~= 0);
    len = sum(~isnan(th(i, :)));
    slope_change_freq(i) = b / len;
end
    
figure('pos', [50 50 100 100]); 
h = histogram(slope_change_freq * 30 * 10, 'binwidth', 0.5, 'normalization', 'probability');
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 1; h.LineWidth = 1.5;
xlim([-0.3 8])
xlabel('N_{reverse-turns}')
setfigprop(700, 800, 4, 40)

% Concatenate all the slopes
meanSlope = [];
for c = 1:57
    s = slope{c};
    meanSlope(c) = nanmean(abs(s));
end

%% Histogram the angular velocity (of all the sections)
rect = [100 100 600 600];
% figure('pos', rect); 
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(meanSlope, 'BinWidth', 0.7, 'normalization', 'probability'); 
% h.FaceColor = 'k'; h.FaceAlpha = 0.7; 
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0.1; h.LineWidth = 2;
xlabel('\omega (rad/h)')
% ylabel('n')
% ylim([0 20])
% title('angular speed')
% set(gca, 'fontsize', 45, 'linewidth', 4.5)
box off
setfigprop(500, 600, 4, 40)

%% Histogram the rotation period (inverse of angular velocity)
rect = [100 100 600 600];
% figure('pos', rect); 
new_fig([50 50 500 600], 4.5, 45, [], [], [], [], "");
h = histogram(2*pi ./ meanSlope * 60, 'BinWidth', 60*1.5, 'normalization', 'probability'); 
% h.FaceColor = 'k'; h.FaceAlpha = 0.7; h.EdgeAlpha = 0;
h.FaceColor = 0.05*[1 1 1]; h.EdgeAlpha = 0.1; h.LineWidth = 1.5;
xlabel('{\itT} (min)')
% ylabel('n')
% xlim([0 12])
% ylim([0 20])
xlim([0 810])
xt = 0:240:900; xticks(xt)
scatter(nanmean(2*pi ./ meanSlope * 60), 0.0, 350, [1 0 0], '^', 'filled')
% title('angular speed')
% set(gca, 'fontsize', 35, 'linewidth', 3.5)
box off
setfigprop(500, 440, 4, 40)

%% Step 5. Now, I can 'flatten' the theta vs. time
% Make all negative slopes from 'Mnew' positive, and cumulative sum to get
% the angles.
thFlat = nan(size(th));
time = nan(size(th));
for c=1:57
    s = diff(th(c,Mnew{c}));
    dtheta = [];
    ssign = sign(s(1));
    for si = 1:length(s)
        if sign(s(si)) ~= ssign 
            w = - dth(c, Mnew{c}(si) : Mnew{c}(si+1));
        else
            w = dth(c, Mnew{c}(si) : Mnew{c}(si+1));
        end
        dtheta = [dtheta w];
    end
 
%     as = ssign * abs(s);
    theta = cumsum(dtheta);
    theta = [0 theta];
    thFlat(c, 1:length(theta)) = theta;
    t = (Mnew{c}-1) / 30; % hour
    time(c, 1:length(t)) = t;
end


% Display the result
T = size(th, 2);
figure; hold on 
% sample = 43;
sample = 27;

for c=1:58
    subplot(7,9,c)
    plot((0:T-1)/30, thFlat(c,:), 's-', 'linewidth', 1.5); hold on
    plot((0:T-1)/30, th(c,:), 'linewidth', 1.5)
    for m = Mnew{c}
        plot((m-1)/30, th(c,m), 'x', 'linewidth', 2)
    end
    ylim([-40 40])
    title(num2str(c))
end

%% All new flattened angles on a graph
figure('pos', [100,100,694,600]); hold on 
for c=1:57
    plot((0:T-1)/30, thFlat(c,:), 'linewidth', 2.4, 'color', [rand(1,3) 0.6]);
end
r = refline(0, 0); r.Color = 'k'; r.LineWidth = 3; r.LineStyle = '--';
xlim([0 20])
xlabel('t (min)')
ylabel('\theta (rad)')
% ylim([-60 60])
ylim([-16*pi 16*pi])
set(gca, 'fontsize', 37, 'linewidth', 3.7)
yt = [-16*pi : 8*pi : 16*pi]; yticks(yt);
ytl = {'-16\pi' '-8\pi' '0' '8\pi' '16\pi'}; yticklabels(ytl)
xt = 0:5:20; xticks(xt)
xt = xt * 60; xticklabels(xt)
box off
setfigprop(1200, 880, 4, 40)

% Just one sample

%% 3. Display the dth(left), D(right) vs time : PAIRS
p = [29, 43, 33];
p = 31;
T = size(dth, 2);
i = 1;
len = 20;
% col = {'r', 'b', [0 0.7 0], [0.93 0.69 0.13]};
fig = figure('pos', [100 100 900 400]); 
for pind = p
%     subplot(3,1,i)
%     subplot(1, 1, i)
%     f = abs(dth(pi, :)) .* D(pi,1:end-1) * sc;  % Tangential speed
    f = abs(dth(pind, :)); % Only the 'w'
%     plot(0:T-2, squeeze(abs(dth(p, n, :)))); ylabel('|\omega|D')
    f_new = nan(T-len,1);
    for t=1:T-len-1
        f_new(t) = mean(f(t:t+len));
    end
    
    yyaxis left
%     plot((0:T-len-1)/(T-len-1)*(T-1)*2/60, squeeze(f_new) * 30, 'linewidth', 4); ylabel('v_{\theta}(\mum/h)')
    plot((0:T-len-1)/(T-len-1)*(T-1)*2/60, squeeze(f_new) * 30, 'linewidth', 4); ylabel('\omega (rad/h)')
%     ylim([0 4.5])
%     ylim([0 60])
%     xlim([2 14])
%     xticklabels(string(0:2:12))
%     yticks(0:20:60)
    yticks(0:2:6)
    
    yyaxis right
    plot((0:T)*2/60, squeeze(D(pind ,:)) , 'linewidth', 4); ylabel('D(\mum)')
%     ylim([5 30])
%     xlim([2 14])
%     plot((0:T)*2/60, shape, 'linewidth', 4); ylabel('shape index')
    
%     xticklabels(string(0:2:12))
    
    i = i + 1;
%     set(gca, 'fontsize', 25, 'linewidth', 3)
    box off
end

% xlim([0 18.4])
xlim([5 15])
set(gca, 'fontsize', 35, 'linewidth', 3.5)
xlabel('hour')

%% 4. Pearson correlation (histogram)
% X = D(:, 1:end-1);
% Y = abs(dth(:, :)) * 30; % Rad per hour

X = newD(:, 1:end-1);
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

%% Pearson coefficient histogram 
figure('pos', rect); 
h = histogram(r, 'binwidth', 0.075);
h.FaceColor = 'k'; h.FaceAlpha = 0.7; h.EdgeAlpha = 0;
% xlim([-1 1])
xlabel('\rho')
ylabel('n')
xlim([-0.6 0.1])
% title('Quad : time window = full')
% title('time window = full')
% title('Pearson correlation coefficient')
set(gca, 'fontsize', 35, 'linewidth' ,3.5)
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5. "Flattening" the angles-vs-time 
% Flattens virtually 'curved' angles so that I can obtain 'pure slope'
T = size(th,2)-1;
figure; 
plot(0:T, th, 'linewidth', 1.5)
figure; 
histogram(dth(:))
%% Angular speed histogram (each pair)
% cftool

% edges = -0.5 : 0.03 : 0.5;
edges = 0 : 0.001 : 0.5;
wHist = nan(length(edges)-1, 57);
for c = 1:57
    cnt = histcounts(abs(dth(c, :)), edges );
    wHist(:, c) = cnt;
end

meanwHist = nanmean(wHist, 2);
stdHist = nanstd(wHist, 0, 2);


figure; shadedErrorBar(edges(1:end-1), meanwHist , stdHist)
hold on ; plot(edges(1:end-1), meanwHist, '.-', 'linewidth', 5, 'markersize', 40)
ylim([-1 8])
xlim([0 0.1])
xlabel('|w|')

% figure; plot(edges(1:end-1), wHist, '.-', 'markersize', 10, 'linewidth', 1)


%% Smooth the curves (smoothing splines)
% figure;
wb = waitbar(0);
interval = 1/10;
newTh = nan(size(th,1), size(th,2)/interval);
for c = 1:57
    y = th(c, :)';
    x = 1:668; x = x';
    notNaN = ~isnan(y);
    y = y(notNaN);
    x = x(notNaN);
    f = fit(x, y, 'smoothingspline', 'smoothingparam', 1e-4);
    
    % Slice the x in to more pieces
    newx = 1:interval:x(end);
    newTh(c, 1:length(newx)) = f(newx);
%     subplot(7,9,c)
%     plot(f, x, y)
%     legend('off')
%     pause(0.5)
    waitbar(c/57)
end
close(wb)
%% newTh's angles
newT = 1:length(newTh);

% figure; plot(newT, newTh)
figure;
for c = 1:57
    subplot(7,9,c)
    plot(newT, newTh(c, :), 'linewidth', 2)
    title(num2str(c))
    refline(0, 0)
    ylim([-40 40])
end

%% new angular speed
newW = diff(newTh, 1, 2);
newDw = diff(newTh, 2, 2);

% figure; histogram(newW(:))

edges = -0.03:0.002:0.03;
% edges = 0 : 0.002 : 0.03;
edgesx = edges(1:end-1) + (edges(2) - edges(1))/2;
newHist = nan(c, length(edges)-1);
meanEdge = nan(c, 1);
stdEdge = nan(c, 1);
cvEdge = nan(c, 1);

for c = 1:57
    h = histcounts((newW(c, :)), edges, 'normalization', 'probability');
    Total = sum(h);
    Mean = sum(h.*edgesx)/Total;
    meanEdge(c) = Mean;
    Std = sqrt(sum(((edgesx - Mean).^2).*h)/Total);
    stdEdge(c) = Std;
    cvEdge(c) = Std / Mean;
    newHist(c, :) = h;
end


%%
figure;
for c=1:16
    subplot(4,4,c)
%     plot(edgesx, newHist(c,:), 'linewidth', 2)

    

%     yyaxis left
%     plot(newW(c, :))
%     ylim([-0.04 0.04])
%     
% %     title(num2str(meanEdge(c)))
% %     title(num2str(stdEdge(c)))
%     
%     cv = cvEdge(c);
%     cv = round(cv, 3);
% %     ti = title(num2str(cv));
%     title(num2str(c))
%     if cv < 0.5
%         ti.Color = 'r';
%     end
% %     l = line([0 0], [0 0.4]);
% %     l.Color = 'r';
% %     l.LineWidth = 1;
% %     l.LineStyle = ':';
%     
%     refline(0, 0)
%     
%     
%     yyaxis right
%     plot(newTh(c, :))
%     ylim([-40 40])


    yyaxis left
    plot(newTh(c, :), 'linewidth', 2)
    ylim([-40 40])
    refline(0, 0)
    
    yyaxis right
    plot(newDw(c, :))
    ylim([-5e-4 5e-4])
    title(num2str(c))
end

%%



%%
figure; plot(newTh(1, :))


figure; plot(edgesx, newHist(:, :))



%%
thres = 0.1;
dth2 = dth;
% B = nanmean(newDth, 3);
% thr = abs(B) * 1;
% newDth(dth < A*thr) = abs(dth((dth < A*thr)));
A = 2*(nanmean(dth2, 2) > 0)-1;
Ap = A > 0;
An = A < 0;
I1 = logical((dth < -thres).*Ap);
I2 = logical((dth > thres).*An);
dth2(I1) = - dth(I1);
dth2(I2) = - dth(I2) ;

sz = size(dth2);
sz(2) = sz(2)+1;
th2 = nan(sz); % all slopes are positive 
th2(:,1) = 0;
for i=1:size(position,2)
    a = dth2(i,:);
    a = squeeze(a);
    a = a(~isnan(a));
    ac = cumsum(a);
    th2(i,2:length(ac)+1) = ac; % Start time doesn't matter much
end



%% Check
figure; hold on 
T = size(th, 2) - 1;
for i=1:size(position,2)
%     subplot(6, 10, i)
    subplot(7, 9, i)
    for j=1
        plot((0:T)*2/60, squeeze(th2(i,:)), 'linewidth', 2); hold on
        plot((0:T)*2/60, squeeze(th(i,:)), 'linewidth', 2, 'color', 'r');
%         plot(squeeze(th3(j,i,:)), 'x', 'linewidth', 2, 'color', 'r');
        xlim([0 15])
        ylim([-30 30])

        ti = title(num2str(i));
%         if ismember(i, noRedir)
%             ti.Color = [1 0 0];
%         elseif ismember(i, oneRedir)
%             ti.Color = [0 1 0];
%         elseif ismember(i, wobble)
%             ti.Color = [0 0 1];
%         end

    end
    r = refline(0, 0); r.LineStyle = ':'; r.Color = [0 0 0];
    refline(0, pind)
    refline(0, -pind)
    refline(0, 2*pind)
    refline(0, -2*pind)
end
xlabel('hours')
ylabel('rad')


%%

% Plot all flattened angles
figure; 
plot((0:T-1)*2/60, th2, 'linewidth', 1.5)
xlim([0 12])
set(gca, 'fontsize', 25, 'linewidth', 2.5)
xlabel('hour')
ylabel('rad')
xticks(0:6:12);
box off
yt = -12*pind:6*pind:12*pind;
yticks(yt)
% ytl = {'-12\pi','-8\pi', '-4\pi', '0', '4\pi', '8\pi', '12\pi'};
ytl = {'-12\pi', '-6\pi', '0', '6\pi', '12\pi'};
yticklabels(ytl)



% Fit the theta to linear functions
A = 1:57;
a = [10 14 16 54]; 
A = setdiff(A, a);
th2(a, :) = th(a, :); % Use these cells' original theta

figure;
slopes = nan(length(A), 2);
for i=1:length(A)
    y = squeeze(th2(A(i), :));
    y = y(~isnan(y));
    y = abs(y);
    x = (0:(length(y)-1))*2/60; % In hour
    f = fit(x', y', 'poly1');
    slopes(i,:) = coeffvalues(f);
    plot(f, x, y)
    pause(0.2)
end

figure; histogram(slopes(:,1), 'binwidth', 0.7) 

%% mean "pair correlation"
v = diff(pos_pair, 1, 4);
av = absvec(v , 3);
av(av==0) = j;
u = v ./ av;
u = real(u);
u = squeeze(u);

pc = dotp(u(1,:,:,:), u(2,:,:,:), 3);
mpc = nanmean(pc, 2);
mpc = nanmean(mpc, 4);





























