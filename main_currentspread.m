
close  all
clear all


% electrode parameters
%p.x1 = -2600; p.x2 = 2600; % separation of the electodes, max value
p.x1 = -1000; p.x2 = 1000;
p.lim = [abs(p.x1)*1.5, 500, 1100]; % how much retina to simulate
ss = 25; % simulation resolution
p.y = 0; p.rad = 225/2;
p.gamma = 1.69; p.amp_min =  50;

% plot params
p.plot = 1;
p.pos = [ 42   1088  939   231];

% how bright the current is on the graphs
p.Isc = .4;

rd_v = -p.lim(1):ss:p.lim(1);
z_v = -p.lim(2):ss:p.lim(2);
zv = -p.lim(3):ss:p.lim(3);
[p.X,p.Y,p.Z] = meshgrid(rd_v,z_v,zv);
p.threshOrig = 23.392 *8.5;

k_vals = [3:.25:7];% 3::7; % controls current spread
rd_vals = 0:50:1000; % threshold elevation due to retinal damage
z_vals = 0:50:1000; % different heights

amp_no_rd = NaN(length(k_vals), length(rd_vals));

%% find thresholds as a function of height, where no retinal damage, single electrode!
%%single electrode
for k = 1:length(k_vals)
    p.k = k_vals(k);
    for z = 1:length(z_vals)
        p.z = -z_vals(z);
        p.R1 = sqrt(((p.X-p.x1).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
        p.R1 = p.R1-p.rad; p.R1(p.R1<0) = 0; p.R1 = p.R1./(max(p.R1(:)));
        p.R2 = sqrt(((p.X-p.x2).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
        p.R2 = p.R2-p.rad; p.R2(p.R2<0) = 0; p.R2 = p.R2./(max(p.R2(:)));
        
        p.a1 = 400; p.a2 = 0; 
        p.thresh = p.amp_min; % each electrode goes up to 127
        % find the threshold for the z
        [p,~] = fit('fit_currentspread',p, {'a1'});
        amp_no_rd_single(k, z) = p.a1;
    end
    figure(100); 
    y = amp_no_rd_single(k, :);
    plot(z_vals, y , 'k'); hold on
    xlabel('z'); ylabel('threshold non damaged retina  ')
    text(z_vals(end-1), y(end-1), ['k=',num2str(k_vals(k))]);
end



%% Panel A
p.k = 5; p.z = -500; rd = 200;
p.R1 = sqrt(((p.X-p.x1).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
p.R1 = p.R1-p.rad; p.R1(p.R1<0) = 0; p.R1 = p.R1./(max(p.R1(:)));
p.R2 = sqrt(((p.X-p.x2).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
p.R2 = p.R2-p.rad; p.R2(p.R2<0) = 0; p.R2 = p.R2./(max(p.R2(:)));
p.thresh = 400; 
p.a = 400;p = rmfield(p, 'a1'); p = rmfield(p,'a2');

[p,err] = fit('fit_currentspread',p, {'a'});
% now double it and calculate what the current spread looks like
p.a = p.a*2;
[I, p, Vq] = create_currentspread(p);

figure(1); clf; p.pos = [ 42   1088   939   231];
set(gcf, 'Name', ['k = ', num2str(p.k), ' height = ', num2str(p.z), 't = ', num2str(p.thresh)]);
create_currentspreadfig(I, p); drawnow

        
%% now calculate surfaces for wide range of retinal damage, no figure
thresh = NaN(length(k_vals), length(z_vals), length(r_vals));
amp = thresh; amp_underelectrode = thresh;amp_max = thresh;amp_midpoint = thresh; gv = thresh;
for k = 1:length(k_vals)
    p.k = k_vals(k);
    for z = 1:length(z_vals)
        p.z = -z_vals(z);
        p.R1 = sqrt(((p.X-p.x1).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
        p.R1 = p.R1-p.rad; p.R1(p.R1<0) = 0; p.R1 = p.R1./(max(p.R1(:)));
        p.R2 = sqrt(((p.X-p.x2).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
        p.R2 = p.R2-p.rad; p.R2(p.R2<0) = 0; p.R2 = p.R2./(max(p.R2(:)));
        for r = 1:length(rd_vals)
            p.a = 400;
            p.thresh = rd_vals(r) +  amp_no_rd(k, z);
            thresh(k,z,r) = p.thresh;
            % find the amplitude to reach this threshold
            [p,err] = fit('fit_currentspread',p, {'a'});
            amp(k,z,r) = p.a;
            
            % now double it and calculate what the current spread looks like
            p.a = p.a*2;
            [I, p, Vq] = create_currentspread(p); 
          
            % current value directly under the electrode on retinal surface
            amp_underelectrode(k, z, r) = interpn(unique(p.Y), unique(p.X), unique(p.Z), I, 0, p.x2, 0);
            % max current value on the surface of the retina
            amp_max(k, z, r) = max(max(I(:, :, round(size(I, 3)/2))));
            % current value directly in between the two electrodes
            amp_midpoint(k,z,r) = interpn(unique(p.Y), unique(p.X), unique(p.Z), I, 0, 0, 0);
        end
    end
end

%% Dip, Panel B
dip_range = [30 70];
thr_range = [300 500];

dip = NaN(length(k_vals), length(z_vals));
figure(2); clf

for k = 1:length(k_vals)
    dip(k, :)= 100.*(amp_max(k, :, 1)-amp_midpoint(k, :, 1))./amp_max(k, :, 1);
    dip_intersect = interp1(dip(k, :), z_vals, [dip_range(1) min([max(dip(k,:)) dip_range(2)])]);
    ind = z_vals<=dip_intersect(1) & z_vals>=dip_intersect(2);
    plot(z_vals, dip(k,:), 'k'); hold on
    plot(z_vals(ind), dip(k,ind), 'r'); hold on
    patch([0 max(z_vals)  max(z_vals) 0 ], [dip_range(1) dip_range(1) dip_range(2) dip_range(2)] ,[.8 .8 .8],  'FaceColor',[ .8 .8 .8], 'EdgeColor', 'none', 'FaceAlpha',.3);
end
set(gca, 'YLim', [0 100]);
xlabel('z (electrode-retina distance')
ylabel(' 1 percept  <-   dip size  ->   2 percepts ')

%% Intersection of constraint curves, Panel C
figure(3); clf
for k = 1:length(k_vals)
    dip_intersect = interp1(dip(k, :), z_vals, [dip_range(1) min([max(dip(k,:)) dip_range(2)])]);
    for r = 1:length(rd_vals)
        goodvals(k, :, r)= z_vals<=dip_intersect(1) & z_vals>=dip_intersect(2) & amp(k,:,r)>=thr_range(1) & amp(k,:,r)<=thr_range(2);
        plot(z_vals,  amp(k,:,r), 'k'); hold on
        plot(z_vals(goodvals(k,:,r)),  amp(k, goodvals(k, :, r),r), 'r'); hold on
    end
end
xlabel('z (electrode-retina distance')
ylabel('threshold (uAmps)')
patch([0 max(z_vals)  max(z_vals) 0 ], [thr_range(1) thr_range(1) thr_range(2) thr_range(2)] ,[.8 .8 .8],  'FaceColor',[ .8 .8 .8], 'EdgeColor', 'none', 'FaceAlpha',.3);

%% patch of good values, Panel D
cmap = hot(length(k_vals)+8);
rd_amp  = amp-amp(:,:,1);
z_amp = amp(:,:,1)-min(amp(:));
z_amp = repmat(z_amp,[1 1  size(rd_amp,3)]);

figure(4); clf

for k=length(k_vals):-1:1
    gv =goodvals(k,:,:);gvi = find(gv(:));
    xv = rd_amp(k, :, :);  xv = xv(:); xg = xv(gvi);
    yv = z_amp(k,:, :); yv = yv(:); yg = yv(gvi);

    if length(xg)>3
        DT = delaunayTriangulation(xg,yg);
        C = convexHull(DT);
        h(k)=patch(DT.Points(C,1),DT.Points(C,2), cmap(k,:), 'EdgeColor', 'none', 'FaceAlpha',.3);
    end
end
xlabel('threshold elevation due to retinal degeneration')
ylabel('threshold elevation due to array lift')



