%% Cluster analysis : recursively finding neighbors
% Serially does clustering analysis on multiple independent trajectories under
% different parameters 

%%
rect = [100 100 610 550];
threshR = 17;                                               % Threshold distance (simulation scale)
threshAngle = cos(20 / 180 * pi);       % Threshold angle
BoxSize = 300;                                         % Simulation box width (height)
T = size(pos, 1);
time = 100-1 : 40 : T-1;
Cells = size(pos, 3);

%%
clusterN = nan(length(time), Cells, 121);
scan_i = 1;

angles = 10:5:60;
scanRs = 11:1:20;

for angle = threshAngle
    for scanR = threshR
        for param = 1:121
            load(['pos_' num2str(param) '.mat'])

            % Periodic BC
            v = diff(pos, 1);
            pos_per = pos;          % Position of the cells considering periodic boundary conditions
            wb = waitbar(0, 'Periodic BC');
            for t=1:T-1
                pos_per(t+1,:,:) = pos_per(t, :,:) + v(t,:,:);
                A = pos_per(t+1, :, :);
                pos_per(t+1,:,:) = A -BoxSize*(A > BoxSize) + BoxSize*(A < 0);
                waitbar(t/T)
            end
            close(wb)

            % Get unit vector
            av = sqrt(v(:, 1, :, :).^2 + v(:,2,:,:).^2);
            av(av==0) = 1i;
            u = real(v ./ av);
            clear av v

            % Recursive clustering
            tt = 1; 
            wb = waitbar(0, ['clustering (param=' num2str(param)]);
            for t = time

                remainingIdx = 1:Cells;
                selected = [];

                N = 1;
                while length(selected) < Cells
                    remainingIdx = setdiff(remainingIdx, selected);
                    indexStack = [];
                    indexQueue = [remainingIdx(1)];
                    count = 0;
                    while ~isempty(indexQueue)
                        i = indexQueue(1);
                        indexQueue(indexQueue == i) = [];
                        selected = union(selected, i);
                        indexStack = union(indexStack, i);

                        xy = squeeze(pos_per(t, :, :));
                        unit = squeeze(u(t, :, :));
                        idx = findMemberIdx(xy, unit, i, BoxSize, scanR, threshAngle); % Does not include itself
                        idx  = setdiff(idx, selected); % Remove already selected cells 
                        selected = union(selected, idx);
                        indexQueue = union(indexQueue, idx);
                        indexStack = union(indexStack, idx);

                    end

                    clusterN(tt, N, param) = length(indexStack);


                    N = N + 1;
                end
                tt = tt + 1;
                waitbar(tt / length(time))
            end
            close(wb)

        end

        scan_i = scan_i + 1;

    end
end
%%
%% Verity (test)
% Total number of cells

s = nansum(clusterN, 2);

% Display
figure;
for i = 1:length(clusterIndex)
    idx = clusterIndex{i};
    x = pos(t, 1, idx);
    y = pos(t, 2, idx);
    ux = u(t, 1, idx);
    uy = u(t, 2, idx);
    
    q = quiver(x, y, ux, uy, 'color', [rand(1, 3)], 'linewidth', 2, 'autoscale', ...
        'off', 'markersize', 20, 'maxheadsize', 0.1);
    hold on 
end

%% Show heatmap of the clusterN
mclusterN = nanmean(clusterN, 2);
mclusterN = nanmean(mclusterN, 1);
mclusterN = squeeze(mclusterN);
mclusterN = reshape(mclusterN, 11, 11);

figure('pos', rect); 
h = heatmap(mclusterN(end:-1:1, :), 'gridvisible', 'off');
h.XDisplayLabels = string(S);
h.YDisplayLabels = string(E(end:-1:1));
h.CellLabelColor = 'none';
% h.Colormap = parula;
colormap(parula);
set(gca, 'fontsize', 12)
xlabel('S')
ylabel('E')





