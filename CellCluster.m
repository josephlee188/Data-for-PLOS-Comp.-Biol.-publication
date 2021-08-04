% [Reference] The function "tube_surf" uses MATLAB function "tubeplot" written by Anders Sandberg, asa@nada.kth.se, 2005

%%
classdef CellCluster
    properties (Access = public)
        blu = [0 0.5 0.8];
        ora = [1.00,0.45,0.35];
        yel = [0.93 0.69 0.13];
        pur = [0.60,0.30,0.70];
        gra = [0.2 0.2 0.2];
        gre = [0.40,0.65,0.35];
        sky = [0.30,0.75,0.93];
        sky2 = [0.68 0.85 0.90];
    end
    
    properties (Access = public)
        pos;
        T;
        N;
        n;
        scale = 1000/340.5;
        colors;
    end
    
    methods (Access = public)
        function obj = CellCluster(pos)
            obj.pos = pos;
            obj.T = size(pos, 4);
            obj.N = size(pos, 2);
            obj.n = size(pos,1);
            obj.colors = cat(1, obj.blu, obj.ora, obj.yel, obj.pur, obj.gre, obj.sky, obj.gra, obj.sky2);
        end
        
        function dth = ang_vel(obj)
            mp = nanmean(obj.pos, 1);
            dv = obj.pos - mp;
            adv = absvec(dv, 3);
            adv(adv==0) = 1j;
            du = dv./adv; 
            du = real(du);
            
            dth = nan(obj.n, obj.N, obj.T-1);            
            for t = 1:obj.T-1
                u1 = du(:,:,:,t);
                u2 = du(:,:,:,t+1);
                crossu = u1(:,:,1).*u2(:,:,2) - u1(:,:,2).*u2(:,:,1);
                crossu = - crossu;
                ang = asin(crossu);
                dth(:,:,t) = ang;
            end
        end
        
        function D = cell_dist(obj)
            mp = nanmean(obj.pos, 1);
            D = obj.pos - mp;
            D = absvec(D, 3);
            D = squeeze(D) * obj.scale; % Calibrated
        end
        
        function th = ang(obj)
            dth = ang_vel(obj);
            th = nan(obj.n, obj.N, obj.T);
            th(:,:,1) = 0;
            for ii=1:obj.n
                for i=1:obj.N
                    a = dth(ii, i, :);
                    a = squeeze(a);
                    a = a(~isnan(a));
                    ac = cumsum(a);
                    th(ii, i, 1:length(ac)) = ac; % Start time doesn't matter much
                end
            end
        end
        
        function init_th = initial_th(obj, t)
            mp = nanmean(obj.pos, 1);
            dv = obj.pos - mp;
            adv = absvec(dv, 3);
            adv(adv==0) = 1j;
            du = dv./adv;
            du = real(du);
            
            x_unit = zeros(size(du));
            x_unit(:,:,1,:) = 1;
            cross = du(:, :, 1, t).*x_unit(:, :, 2, t) - du(:, :, 2, t).*x_unit(:, :, 1, t);
            dot = du(:, :, 1, t).*x_unit(:, :, 1, t) + du(:, :, 2, t).*x_unit(:, :, 2, t);
            init_th = sign(cross).*acos(dot);
        end
            
        function result = smoothing(obj, input, sp)
            % input : (number of cells, number of samples, time)
            result = nan(size(input));
            for ii = 1:obj.n
                for i = 1:obj.N
                    y = input(ii, i,:); y = squeeze(y);
                    idx = find(~isnan(y));
                    if length(idx) > 2
                        x = 1:length(y); x = x';
                        x = x(idx);
                        y = y(idx);
                        if ~isempty(y)
                            f = fit(x, y, 'smoothingspline', 'SmoothingParam', sp);
                            result(ii, i, idx) = reshape(f(x), 1, 1, length(idx));
                        end
                    end
                end
            end
        end
        
        function result = subtract_smoothing(obj, input, sp)
            % input : (number of cells, number of samples, time)
            s = smoothing(obj, input, sp);
            result = input - s;
        end

        function t_period = get_period_fft(obj, input)
            intv = 2/60; % Interval (hour)
            Fs = 1/intv;
            t_period = nan(obj.n, obj.N);
            j = 1;
            for ii = 1:obj.n
                for i=1:obj.N
                    dd = squeeze(input(ii, i, :));
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
                        % Smooth the curve
                        f = fit(t(2:end)', P1(2:end), 'smoothingspline', 'SmoothingParam', 0.90);
                        f = f(t(2:end));
                        ind = find(max(f) == f) + 1;
                        t_period(ii, i) = t(ind);
                        j = j + 1;
                    end
                end
            end
        end
        
        function pcorrel = pearson_correlation(obj, input1, input2)
            input1 = input1(:, :, 1:obj.T-1);
            input2 = input2(:, :, 1:obj.T-1);
            m1 = nanmean(input1, 3);
            m2 = nanmean(input2, 3);
            A = (input1 - m1) .* (input2 - m2);
            A = nansum(A, 3);
            B = sqrt(nansum((input1 - m1).^2, 3) .* nansum((input2 - m2).^2, 3));
            pcorrel = A(:) ./ B(:);
        end
        
        function plot3d_traj(obj, cell_ind, timespan, rect, lim, viewangle, col)
%             blu = [0 0.45 0.74];
%             ora = [0.85 0.33 0.1];
%             yel = [0.93 0.69 0.13];
%             pur = [0.49 0.18 0.56];
%             color = cat(1, blu, ora, yel, pur);

            position = squeeze(obj.pos(:, cell_ind, :, timespan));
            if length(size(position))<3; position = reshape(position, [1 size(position)]); end
            mp = nanmean(position, 1);
            mp = squeeze(mp);
%             figure('pos', rect); hold on
            z = timespan - timespan(1);
            z = z / 30;
            for i = 1:obj.n
                x = (position(i,1,:) - mp(1,1))* obj.scale; x = squeeze(x);
                y = (position(i,2,:) - mp(2,1))* obj.scale; y = squeeze(y);
                % Smooth the traj
                indx = ~isnan(x); indy = ~isnan(y);
                ind = indx & indy;
                if sum(ind) > 2
                    x_fit = fit(z(ind)', x(ind), 'smoothingspline', 'smoothingparam', 0.98);
                    y_fit = fit(z(ind)', y(ind), 'smoothingspline', 'smoothingparam', 0.98);
                    x_smooth = x_fit(z(ind));
                    y_smooth = y_fit(z(ind));
                    % Draw 3D traj
                    if ~isempty(col)
                        plot3(x_smooth, y_smooth, z(ind), 'linewidth', 6, 'color', [col 0.9]);
                    else
                        plot3(x_smooth, y_smooth, z(ind), 'linewidth', 6, 'color', [obj.colors(i,:) 0.9]);
                    end
                    % Draw 2D projections
                    plot3(x_smooth(1), y_smooth(1), 0.0, 'd', 'markersize', 20, ...
                        'markerfacecolor', obj.colors(i,:), 'markeredgecolor', obj.colors(i,:))
                    %     plot3(x_smooth, zeros(size(x_smooth))-lim, z(ind), 'linewidth', 6, 'color', ones(1, 3)*0.7);
%                     plot3(zeros(size(y_smooth))-lim + 0.5, y_smooth, z(ind), 'linewidth', 6, ...
%                         'color', [obj.colors(i, :) 0.3]);
                    plot3(x_smooth, zeros(size(x_smooth))-lim + 0.5, z(ind), 'linewidth', 6, ...
                        'color', [obj.colors(i, :) 0.3]);
                    %     if i == 1
                    %         plot3(x_smooth, y_smooth, zeros(size(x_smooth)) + 14, 'linewidth', 6, 'color', ones(1, 3)*0.7);
                    %     end
                end
            end
            grid on
            
            xlim([-lim lim])
            ylim([-lim lim])
            
            set(gca, 'fontsize', 45, 'linewidth', 4.5, 'cameraposition', viewangle);
            xt = -40:20:40;
            yt = xt;
            xticks(xt);
            yticks(yt);
            
%             xlabel('x (\mum)')
%             ylabel('y (\mum)')
            zlabel('t (hour)')
            
        end

        function overlay_image(obj, img, cell_ind, lim, offsetx, offsety)
            position = squeeze(obj.pos(:, cell_ind, :, 1));
            mp = nanmean(position, 1);
            mp = squeeze(mp);
            % Overlay the first frame at the bottom
            w = 25; scaleImg = 1.2;
            x = mp(1); x = round(x);
            y = mp(2); y = round(y);
            x1 = max(1, x-w); x2 = min(size(img, 2), x+w);
            y1 = max(1, y-w); y2 = min(size(img, 1), y+w);
            im_crop = img(y1:y2, x1:x2);
            
            wid = size(im_crop, 1);
            [X,Y] = meshgrid(0:wid-1, 0:wid-1);
            X = X - wid/2; X = X * scaleImg;
            Y = Y - wid/2; Y = Y * scaleImg;
            X = X * obj.scale + offsetx;
            Y = Y * obj.scale + offsety;
            Z = zeros(size(X))+0.01;
            
            [x, y] = meshgrid(-lim:lim, -lim:lim);
            z = zeros(size(x));
            surf(x,y,z, ones(size(x)),'facecolor', 'texturemap', 'edgecolor', 'none'); hold on
            surf(X, Y, Z, im_crop, 'facecolor', 'texturemap', 'edgecolor', 'none');
            colormap(gray)
        end
        
        function u = get_unit_vector(obj)
            v = diff(obj.pos, 1, 4);
            absv = absvec(v, 3);
            v(v == 0) = 1i;
            u = v ./ absv;
            u = real(u);
        end
        
        function pc = pair_correlation(obj)
            u = get_unit_vector(obj);
            c = u(1, :, 1, :) .* u(2, :, 1, :) + u(1, :, 2, :) .* u(2, :, 2, :);
            pc = squeeze(c);
        end
        
        function qc = quad_correlation(obj)
            % Correlation between each NEIGHBOR pairs in quads
            u = get_unit_vector(obj);
            qc = nan(5, obj.N, obj.T-1);
            for t = 1:obj.T-1
                for j = 1:obj.N
                    for i = 1:5
                        if ~isnan(obj.pos(i, j, 1, t))
                            c = u(i, j, 1, t) .* u(:, j, 1, t) + u(i, j, 2, t) .* u(:, j, 2, t);
                            c(i) = nan;
                            qc(i, j, t) = nanmean(c);
%                             d = obj.pos(i, j, :, t) - obj.pos(:, j, :, t);
%                             d = squeeze(d);
%                             d = sqrt(d(:,1).^2 + d(:,2).^2);
%                             [~, idx] = sort(d);
%                             idx1 = idx(2); idx2 = idx(3);
%                             c1 = u(idx1, j, 1, t) * u(i, j, 1, t) + u(idx1, j, 2, t) * u(i, j ,2, t);
%                             c2 = u(idx2, j, 1, t) * u(i, j, 1, t) + u(idx2, j, 2, t) * u(i, j ,2, t);
%                             qc(i, idx1, j, t) = c1;
%                             qc(i, idx2, j, t) = c2;
                        

                        end
                    end
                end
            end
%             mqc = squeeze(nanmean(qc, 2));
        end
                            
        function result = Msd(obj)
            result = zeros(obj.N, obj.T);
            for dt = 0 : obj.T-1
                t0 = 1:obj.T-dt;
                t1 = t0 + dt;
                d = obj.pos(:, :, :, t1) - obj.pos(:, :, :, t0);
                d2 = d(:, :, 1, :).^2 + d(:, :, 2, :).^2;
                result(:, dt+1) = squeeze(nanmean(d2, 4));
            end
            result = result * obj.scale^2;
        end
        
        function result = Auto_correl(obj)
            u = get_unit_vector(obj);
            result = nan(obj.N, obj.T-1);
            for dt = 0 : obj.T-2
                t0 = 1:obj.T-1-dt;
                t1 = t0 + dt;
                c = u(1, :, 1, t0) .* u(1, :, 1, t1) + u(1, :, 2, t0) .* u(1, :, 2, t1);
                result(:, dt+1) = squeeze(nanmean(c, 4));
            end
        end

        function autoc = Auto_correl_cme(obj)
            autoc = nan(obj.N, obj.T);
            for i = 1:obj.N
                position = squeeze(obj.pos(1, i, :, :));
                position = permute(position, [2,1]);
                [x, ~] = cme(position);
                autoc(i, :) = x;
            end
        end
        
        function ptime = cme_time(obj)
            wb = waitbar(0, 'ptime');
            ctime = nan(obj.N, 2);
            for i = 1:obj.N
                position = squeeze(obj.pos(1, i, :, :));
                position = permute(position, [2,1]);
                [~, c] = cme(position);
                waitbar(i/obj.N)
                ctime(i, :) = c;
            end
            ptime = 1 ./ ctime(:,2) / 30; % in hour
            close(wb)
        end
        
        function plength = cme_len(obj)
            wb = waitbar(0, 'plength');
            clen = nan(obj.N, 4);
            for i = 1:obj.N
                position = squeeze(obj.pos(1, i, :, :));
                position = permute(position, [2,1]);
                [~, c] = cme2(position);
                waitbar(i/obj.N)
                clen(i, :) = c;
            end
            plength = 1 ./ clen(:, [2 4]) * obj.scale;
            close(wb)
        end
        
        function plot_scatter(obj, ms, lw, malpha, color, x, y, ystd)
            eb = shadedErrorBar(x, y, ystd);
            eb.mainLine.LineStyle = 'none';
            eb.patch.FaceColor = color; 
            eb.edge(1).LineStyle = 'none'; 
            eb.edge(2).LineStyle = 'none';
%             plot(x, y, 'linewidth', lw, 'color', color);
            scatter(x, y, ms, ...
                'markerfacecolor', color, 'markerfacealpha', malpha, 'markeredgealpha', 0);
            
        end
        
        function [Mccorrel, Ncount] = correl_map(obj, edges, timespan)
            
            u = get_unit_vector(obj);
            u = permute(squeeze(u), [3, 2, 1]);
            u = u(timespan, :, :);
            po = obj.pos;
            po = permute(squeeze(po), [3,2,1]);
            po = po(timespan, :, :);
            
            ccorrel = nan(length(timespan),obj.N,obj.N);
            distx = nan(size(ccorrel));
            disty = nan(size(ccorrel));
            for i=1:obj.N
                ccorrel(:,i,:) = u(:,1,i).*u(:,1,:)+u(:,2,i).*u(:,2,:);
                dr = po-po(:,:,i);
                disty(:,i,:) = dr(:,1,:).*u(:,1,i) + dr(:,2,:).*u(:,2,i);
                distx(:,i,:) = dr(:,1,:).*u(:,2,i) - dr(:,2,:).*u(:,1,i);
            end
            
%             distx = distx * obj.scale;
%             disty = disty * obj.scale;
%             map = parula(100);
            I = ~isnan(ccorrel(:)) & ~isnan(distx(:)) & ~isnan(disty(:));
            % edges = -700:15:700;
%             edges = -41 : 2 : 41;
            % Nperm = 750; % Maximum number of datasets which comes into a small box
            % Nperm = 2000;
            X = discretize(distx(I),edges);
            Y = discretize(disty(I),edges);
            C = ccorrel(I);
            clear distx disty
            clear ccorrel
            Mccorrel = nan(length(edges),length(edges));
            Ncount = nan(length(edges), length(edges)); % Probability density
            % STDccorrel = nan(length(edges),length(edges));
            ii = 0;
            wb = waitbar(0);
            for y=1:length(edges)
                for x=1:length(edges)
%                     cc = C((X==x)&(Y==y));
%                     cc = cc(~isnan(cc));
                    %         if length(cc) >= Nperm
                    %             R = randperm(length(cc), Nperm);
                    %             Mccorrel(x,y) = nanmean(cc(R));
                    %             Ncount(x,y) = Nperm;
                    %             STDccorrel(x,y) = nanstd(cc(R));
                    %         else
                    %             Mccorrel(x,y) = nanmean(cc);
                    %             STDccorrel(x,y) = nanstd(cc);
                    %             Ncount(x,y) = length(cc);
                    %         end
                    
                    Mccorrel(x,y) = nanmean(C((X==x)&(Y==y)));
                    Ncount(x,y) = sum(~isnan(C((X==x)&(Y==y))));
                    
                    waitbar(ii/length(edges)^2)
                    ii = ii + 1;
                end
            end
            close(wb)  
        end
        
        function [Mccorrel, Ncount] = correl_map_periodic(obj, edges, timespan, celln)            
            u = get_unit_vector(obj);
            u = permute(squeeze(u), [3, 2, 1]);
            u = u(timespan, :, celln);
            po = obj.pos;
            po = permute(squeeze(po), [3,2,1]);
            po = po(timespan, :, celln);
            
            v = diff(po, 1, 1);
            for t = 2:length(timespan)
                po(t, :, :) = po(t-1, :, :) + v(t-1, :, :);
                po(t, :, :) = ...
                    po(t, :, :) - ...
                    300 * (po(t, :, :) > 300) + ...
                    300 * (po(t, :, :) < 0);
            end
            
            ccorrel = nan(length(timespan),length(celln),length(celln));
            distx = nan(size(ccorrel));
            disty = nan(size(ccorrel));
            for i=1:length(celln)
                ccorrel(:,i,:) = u(:,1,i).*u(:,1,:)+u(:,2,i).*u(:,2,:);
                dr = po-po(:,:,i);
                dr = dr - ... 
                    (dr > 150)*300 + ...
                    (dr < -150)*300 ;
                disty(:,i,:) = dr(:,1,:).*u(:,1,i) + dr(:,2,:).*u(:,2,i);
                distx(:,i,:) = dr(:,1,:).*u(:,2,i) - dr(:,2,:).*u(:,1,i);
            end
            
%             distx = distx * obj.scale;
%             disty = disty * obj.scale;
%             map = parula(100);
            I = ~isnan(ccorrel(:)) & ~isnan(distx(:)) & ~isnan(disty(:));
            % edges = -700:15:700;
%             edges = -41 : 2 : 41;
            % Nperm = 750; % Maximum number of datasets which comes into a small box
            % Nperm = 2000;
            X = discretize(distx(I),edges);
            Y = discretize(disty(I),edges);
            C = ccorrel(I);
            clear distx disty
            clear ccorrel
            Mccorrel = nan(length(edges),length(edges));
            Ncount = nan(length(edges), length(edges)); % Probability density
            % STDccorrel = nan(length(edges),length(edges));
            ii = 0;
%             wb = waitbar(0);
            for y=1:length(edges)
                for x=1:length(edges)
%                     cc = C((X==x)&(Y==y));
%                     cc = cc(~isnan(cc));
                    %         if length(cc) >= Nperm
                    %             R = randperm(length(cc), Nperm);
                    %             Mccorrel(x,y) = nanmean(cc(R));
                    %             Ncount(x,y) = Nperm;
                    %             STDccorrel(x,y) = nanstd(cc(R));
                    %         else
                    %             Mccorrel(x,y) = nanmean(cc);
                    %             STDccorrel(x,y) = nanstd(cc);
                    %             Ncount(x,y) = length(cc);
                    %         end
                    
                    Mccorrel(x,y) = nanmean(C((X==x)&(Y==y)));
                    Ncount(x,y) = sum(~isnan(C((X==x)&(Y==y))));
                    
%                     waitbar(ii/length(edges)^2)
                    ii = ii + 1;
                end
            end
%             close(wb)  
        end
            
        function tube_surf(obj, cell_ind, timespan, rect, lim, viewangle, twfactor, interpN, sp, zfactor)
            position = squeeze(obj.pos(:, cell_ind, :, timespan));
            if length(size(position))<3; position = reshape(position, [1 size(position)]); end
            mp = nanmean(position, 1);
            mp = squeeze(mp);
%             figure('pos', rect); hold on
            z = timespan - timespan(1);
%             z = z / 30;
            for i = 1:obj.n
                x = (position(i,1,:) - mp(1,1))* obj.scale; x = squeeze(x);
                y = (position(i,2,:) - mp(2,1))* obj.scale; y = squeeze(y);
                % Smooth the traj
                indx = ~isnan(x); indy = ~isnan(y);
                ind = indx & indy;
                if sum(ind) > 2
                    x_fit = fit(z(ind)', x(ind), 'smoothingspline', 'smoothingparam', sp);
                    y_fit = fit(z(ind)', y(ind), 'smoothingspline', 'smoothingparam', sp);
                    z_new = z(ind); 
                    z_new = z_new - z_new(1);
                    dz = z_new(2) - z_new(1);
                    z_new = z_new(1): dz/interpN : z_new(end); z_new = z_new';
                    x_smooth = x_fit(z_new);
                    y_smooth = y_fit(z_new);
                    z_new = z_new * zfactor;
                    % Draw 3D Tubes
                    vx = diff(x_smooth,1);
                    vy = diff(y_smooth,1);
                    v = sqrt(vx.^2 + vy.^2);
                    tw = v * twfactor; tw = [tw; tw(end)]; tw = [tw tw];
                    tubeplot(x_smooth, y_smooth, z_new, tw, zeros(length(z_new),1), 30);
                    
                    
                    % Draw 2D projections
                    plot3(x_smooth(1), y_smooth(1), 0.0, 'd', 'markersize', 20, ...
                        'markerfacecolor', obj.colors(i,:), 'markeredgecolor', obj.colors(i,:))
                    %     plot3(x_smooth, zeros(size(x_smooth))-lim, z(ind), 'linewidth', 6, 'color', ones(1, 3)*0.7);
%                     plot3(zeros(size(y_smooth))-lim + 0.5, y_smooth, z(ind), 'linewidth', 6, ...
%                         'color', [obj.colors(i, :) 0.3]);
                    plot3(x_smooth, zeros(size(x_smooth))-lim + 0.5, z_new, 'linewidth', 10, ...
                        'color', [obj.colors(i, :) 0.5]);
                    %     if i == 1
                    %         plot3(x_smooth, y_smooth, zeros(size(x_smooth)) + 14, 'linewidth', 6, 'color', ones(1, 3)*0.7);
                    %     end
                end
            end
            
            camlight left
            camlight right
            
            s = findobj(gcf, 'type', 'surface');
            if obj.n == 5
                objn = 4;
            else
                objn = obj.n;
            end
            for i = 1:objn
                ii = objn+1-i;
                s(i).EdgeColor = 'none';
                s(i).FaceColor = obj.colors(ii,:);
                s(i).FaceAlpha = 0.8;
            end
            grid on
            
            xlim([-lim lim])
            ylim([-lim lim])
            
            set(gca, 'fontsize', 45, 'linewidth', 4.5, 'cameraposition', viewangle);
            xt = -40:20:40;
            yt = xt;
            xticks(xt);
            yticks(yt);
            
%             xlabel('x (\mum)')
%             ylabel('y (\mum)')
            zlabel('t (min)')
            dz = z_new(2) - z_new(1); z_new = z_new(1):dz:(z_new(end)+dz*500);
            zt = z_new(1:150*interpN:end);
            z_hour = z_new / zfactor / 30;
            ztl = string(60*z_hour(1:150*interpN:end));
            zticks(floor(zt));
            zticklabels(ztl);
            
        end
        
        function img = remove_boundary(obj, img, bg_black)
            factor = 1.5;
            
            blue2 = uint8(min(obj.blu * 255 * factor, 255));
            orange2 = uint8(min(obj.ora * 255 * factor, 255));
            yellow2 = uint8(min(obj.yel * 255 * factor, 255));
            purple2 = uint8(min(obj.pur * 255 * factor, 255));
            
            white = [255; 255; 255];
            black = [0; 0; 0];
            yellow = [255; 255; 0];        % 1
            blue = [0; 0; 255];                 % 3
            green = [0; 100; 0];              % 2
            orange = [255; 165; 0];       % 4
            
            for i = 2:size(img, 1)-1
                for j = 2:size(img, 2)-1
                    a = squeeze(img(i, j, :));
                    if ~isequal(a, white) && ~isequal(a, black)
                        if isequal(a, yellow)
                            img(i, j, :) = min(obj.blu * 255 * factor, 255);
                        elseif isequal(a, green)
                            img(i, j, :) = min(obj.ora * 255 * factor, 255);
                        elseif isequal(a, blue)
                            img(i, j, :) = min(obj.yel * 255 * factor, 255);
                        elseif isequal(a, orange)
                            img(i, j, :) = min(obj.pur * 255 * factor, 255);
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
                    if ~isequal(a, white') && ~isequal(a, blue2) && ~isequal(a, orange2) && ~isequal(a, yellow2) && ~isequal(a,purple2)
                        b = img(i-1:i+1, j-1:j+1, :);
                        b = reshape(b, 9, 3);
                        b(5,:) = [];
                        
                        wt = sum(b == white', 2);
                        bl = sum(b == blue2, 2);
                        or = sum(b == orange2, 2);
                        ye = sum(b == yellow2, 2);
                        pu = sum(b == purple2, 2);
                        
                        wt_n = sum(wt == 3);
                        bl_n = sum(bl == 3);
                        or_n = sum(or == 3);
                        ye_n = sum(ye == 3);
                        pu_n = sum(pu == 3);
                        
                        [~, ind] = max([wt_n, bl_n, or_n, ye_n, pu_n]);
                        
                        if ind == 1
                            img(i,j,:) = white;
                        elseif ind == 2
                            img(i,j,:) = blue2;
                        elseif ind == 3
                            img(i,j,:) = orange2;
                        elseif ind == 4
                            img(i,j,:) = yellow2;
                        elseif ind == 5
                            img(i,j,:) = purple2;
                        end
                    end
                end
            end
            
            if bg_black == 1
                for i = 1:size(img,1)
                    for j = 1:size(img,2)
                        a = squeeze(img(i, j, :)); a = a';
                        if isequal(a, white')
                            img(i,j,:) = black;
                        end
                    end
                end
            end
            
        end

                
            
        
        
            
                
               
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end

end



% Defining local functions
% functions defined here, outside the class will be visible only inside
% this file





