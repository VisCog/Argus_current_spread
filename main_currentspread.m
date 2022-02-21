
close  all
clear all

lowres = 1;
% electrode parameters
ret.z_range = 0:25:1000; % different heights (microns)


ret.rad = 225/2;
ret.lim = [500, 1500, 1100]; % how much retina to simulate. x, y, z
ret.ss = 25; % simulation resolution

% plot params
plt.pos = [ 42   1088  939   231]; 

% how bright the current is on the graphs
ret.Isc = .4;

xv = -ret.lim(1):ret.ss:ret.lim(1);
yv = -ret.lim(2):ret.ss:ret.lim(2);
zv = -ret.lim(3):ret.ss:ret.lim(3);
[ret.X,ret.Y,ret.Z] = meshgrid(yv,xv,zv); % simulates the retina as a grid

% neurophysiological paramers

ret.k_range = 6:3:20; % controls current spread
ret.a_range = 1:.5:3;
good_ak = [0     1     1     1     1
     1     1     1     1     1
     1     1     1     1     1
     1     1     1     0     0
     1     1     0     0     0]

ret.t_ret_min = 50; % minimum retinal current required to reach perceptual

% threshold
ret.t_rd_range = 0:100:600; % threshold elevation due to retinal damage (mAmps)

fitParams.tol = 0.001;
fitParams.lo = 0;
fitParams.hi = 5000;
fitParams.thr = 1;
fitParams.nreps = 20;

%% calculate thresholds as a function of height, where no retinal damage, single electrode
ret_s = ret; ret_s.t_ret = ret.t_ret_min;
ret_s.x = 0; ret_s.y = 0;
ret_s.t = 400; % current on the electrode

for z = 1:length(ret_s.z_range)
    disp(['fitting z = ', num2str(z),  ' out of ', num2str(length(ret_s.z_range))]);
    ret_s.z = -ret_s.z_range(z);
    ret_s = cs.calc_dist_from_electrode(ret_s);
    for a = 1:length(ret_s.a_range)
        disp(a)
        for k = 1:length(ret_s.k_range)
            ret_s.k = ret_s.k_range(k);
            ret_s.a = ret_s.a_range(a);
            % find the threshold for the z
            ret_s = cs.fit_currentspreadfast(ret_s, fitParams);
            ret_s.t_fit(a, k, z) = ret_s.t;
        end
    end
end


%% plot thresholds as a function of height, where no retinal damage, single electrode

subplot(1,2,1)
for a = 1:length(ret_s.a_range)
    for k = 1:length(ret_s.k_range )
        figure(1); 
        subplot(1,2,1)
        plot(ret_s.z_range, squeeze(ret_s.t_fit(a,k, :)), 'k'); hold on
        xlabel('z'); ylabel('threshold non damaged retina')
        
        text(ret_s.z_range(end-1),  ret_s.t_fit(1,k, end-1), ['k = ',num2str(ret_s.k_range(k))]);
        if k==3 & a==3
              plot(ret_s.z_range, squeeze(ret_s.t_fit(a,k, :)), 'r'); hold on
        end
        set(gca, 'XLim', [0 1000])
        set(gca, 'XTick',[0:200:1000])
        set(gca, 'YLim', [0 700])
        set(gca, 'YTick',[0:100:700])

        subplot(1,2,2)

        plot(log(ret_s.z_range), squeeze(log(ret_s.t_fit(a,k, :))) , 'k'); hold on
        if k==3 & a==3
              plot(log(ret_s.z_range), log(squeeze(ret_s.t_fit(a,k, :))), 'r'); hold on
        end
        xlabel('z'); ylabel('threshold non damaged retina')
        text(log(ret_s.z_range(end-1)),  log(ret_s.t_fit(1,k, end-1)), ['k = ',num2str(ret_s.k_range(k))]);

        set(gca, 'XTick', log([ 100 1000]))
        set(gca, 'XLim', log([50 2000]))
        set(gca, 'XTickLabel',[100 1000])
        logx2raw

        set(gca, 'YTick', log([10  100 1000]))
        set(gca, 'YLim', log([5 2000]))
        set(gca, 'YTickLabel',[10 100 1000])
        logy2raw
      %  val(a, k) = input(' good = 1, bad = 0 ...' );
    end
end


%% calculate surfaces for pair of electrodes, wide range of retinal damage
ret_e = ret; % example electrode pair
ret_e.k = 5; 
ret_e.a = 2;
ret_e.k = 12; 
ret_e.t_lift = 0;
ret_e.x = [-948 948];
ret_e.y = [0 0];
ret_e.t = [400];
ret_e.t_ret = ret.t_ret_min;

ret_p = ret_e; % save these parameters for the full simulation

%% plot, current spreads for 2 electrodes at different heights

% example electrode
z_list = [-150 -750]
for zz = 1:2
    ret_e.z = z_list(zz);
    figure(zz+1);clf
    ret_e = cs.calc_dist_from_electrode(ret_e);
    [ret_e,~] = fit('cs.fit_currentspread',ret_e, {'t'});
    % now double it and calculate what the current spread looks like
    ret_e.t =  ret.t_ret_min * 2;
    ret_e = cs.create_currentspread(ret_e);
    set(gcf, 'Position', plt.pos);
    set(gcf, 'Name', ['height = ', num2str(ret_e.z)]);
    Sxy = ret_e.I(:, :, round(size(ret_e.I, 3)/2));
    ret_e.amp_max = max(Sxy(:));
    [y, x] = find(Sxy== ret_e.amp_max);
    ret_e.amp_mid = interpn(unique(ret_e.Y), unique(ret_e.X),Sxy, 0, 0);
    ret_e.dip= 100.*(ret_e.amp_max-ret_e.amp_mid)./ret_e.amp_max;
    ret_e.loc = [x y];
    cs.create_currentspreadfig(ret_e); drawnow
end



%% calculate for a range of paramter values
ret_p.t_ret = ret_p.t_ret_min;
for z = 1:length(ret_p.z_range)
    disp(['fitting z = ', num2str(z),  ' out of ', num2str(length(ret_p.z_range))]);
    ret_p.z = -ret_p.z_range(z);
    ret_p = cs.calc_dist_from_electrode(ret_p);
    for k = 1:length(ret_p.k_range)
        for a = 1:length(ret_p.a_range)
            if good_ak(a,k)==1
            disp(['fitting k = ', num2str(k),  ' out of ', num2str(length(ret_p.k_range))]);
            ret_p.k = ret_p.k_range(k);
            for r = 1:length(ret_p.t_rd_range)
                ret_p.t_ret= ret_p.t_rd_range(r) + ret_p.t_ret_min;  % assuming that on the retina more current is needed

                % find the amplitude to reach this threshold
                ret_p.t = 400;
                ret_p = cs.fit_currentspreadfast(ret_p, fitParams);
                ret_p.t_fit(a, k,z,r) = ret_p.t;
                % now double it and calculate what the current spread looks like
                ret_p.t = ret_p.t*2;
                ret_p = cs.create_currentspread(ret_p);

                % max current value on the surface of the retina
                ret_p.amp_max(a, k,z,r) = max(max(ret_p.I(:, :, round(size(ret_p.I, 3)/2))));
                % current value directly in between the two electrodes
                ret_p.amp_mid(a, k,z,r) = interpn(unique(ret_p.Y), unique(ret_p.X), unique(ret_p.Z), ret_p.I, 0, 0, 0);
            end
            end
        end
    end
end


%% dipsize as a function of z and k, for a fixed value of r

ret_p.dip= 100.*(ret_p.amp_max-ret_p.amp_mid)./ret_p.amp_max; % 100 means big dip
for a = 1:length(ret_p.a_range)
figure(5); subplot(2, 4, a)
imagesc(ret_p.z_range, ret_p.k_range, squeeze(2.55*ret_p.dip(a, :,:, 1))); 
colormap(gray(256)); hold on
set(gca, 'XTick',ret_p.z_range); xlabel('z');
set(gca, 'YTick',ret_p.k_range); ylabel('k');
ylabel('k'); 
title('dipsize');
%c = colorbar('Ticks',linspace(1, 255, 5), 'TickLabels',{'0','25','50','75','100'});
contour(ret_p.z_range, ret_p.k_range, squeeze(ret_p.dip(a, :,:, 1)), [30 70], 'g')
end
return

%% intersection of constraints, which combinations of z and k, for a fixed r, are plausible


for k = 1:length(ret_p.k_range)
    tmin = ret_p.t_fit(k,1,1);
    for z = 1:length(ret_p.z_range)
        for r = 1:length(ret_p.t_rd_range)
            ret_p.t_rd_fit(k,z,r) = ret_p.t_fit(k,z,r)-ret_p.t_fit(k,z,1);
            ret_p.t_z_fit(k,z,r) = ret_p.t_fit(k,z,r)-(ret_p.t_rd_fit(k,z,r)+tmin); % what extra current did the lift require
            round(ret_p.t_fit(k,z,r))
            round(tmin+ret_p.t_rd_fit(k,z,r) + ret_p.t_z_fit(k,z,r) )
        end
    end
end

%% dip doesn't change with retinal damage
dip_range = [30 70];
thr_range = [177 597];
cmap = [0 0 0 ; 0 0 1 ; 1 0 0 ;  1 1 1 ];
r=3;
ret_p.gv_thr = zeros(size(ret_p.amp_max));ret_p.gv_dip = ret_p.gv_thr;

ret_p.gv_thr(find(ret_p.t_fit(:)>thr_range(1) & ret_p.t_fit(:)<thr_range(2)))=2;
ret_p.gv_dip(find(ret_p.dip(:)>dip_range(1) & ret_p.dip(:)<dip_range(2)))=1;
ret_p.gv = (ret_p.gv_thr+ret_p.gv_dip)==3;

figure(5); clf
image(ret_p.z_range, ret_p.k_range, 1+ret_p.gv_dip(:,:, r)+ret_p.gv_thr(:,:, r)); colormap(cmap); hold on
set(gca, 'XTick',ret_p.z_range); xlabel('z');
set(gca, 'YTick',ret_p.k_range); ylabel('k');
ylabel('k'); 
title(['intersection ret damage  = ', num2str(ret_p.t_rd_range(r))]);
c = colorbar('Ticks',[1 85 169 255], 'TickLabels',{'neither','dip','thr','dip+thr'});


%% patch of good values, Panel D
cmap = hot(length(ret_p.k_range)+8);
figure(6); clf

for k=length(ret_p.k_range):-1:1
    gv =ret_p.gv(k,:,:); gv_idx = find(gv(:)); % find the values that are ok
    t_z = ret_p.t_z_fit(k, :, :); t_z = t_z(:); t_z = t_z(gv_idx); % collect the increase in thr due to lift
    t_all = ret_p.t_fit(k, :, :);t_all = t_all(:); t_all = t_all(gv_idx);
    t_rd = ret_p.t_rd_fit(k, :, :); t_rd = t_rd(:); t_rd = t_rd(gv_idx); 
    zv = repmat(ret_p.z_range, 1, length(ret_p.t_rd_range)); zv = zv(:); zv = zv(gv_idx);
    round(min(t_z(:)+t_rd(:)))
    round(max(t_z(:)+t_rd(:)))
    round(max(t_all(:)))

    if length(gv_idx)>6
        DT = delaunayTriangulation(t_z,t_rd);
        C = convexHull(DT);
        h(k)=patch(DT.Points(C,1),DT.Points(C,2), cmap(k,:), 'EdgeColor', 'none', 'FaceAlpha',.3); hold on
    end
end
set(gca, 'XLim', [ 0 300])
set(gca, 'YLim', [ 0 600])
xlabel('threshold elevation due to retinal degeneration')
ylabel('threshold elevation due to array lift')


%% thresholds
tdata = {1	'B10'	89.0; 1	'B10'	93.0; 2	'A8'	153; 1	'B9'	177; ...
1	'F10'	177; 2	'A4'	194; 1	'F10'	202; 1	'C6'	210; 2	'F2'	226; ...
1	'E9'	242; 1	'A8'	250; 2	'D1'	274; 2	'F2'	274; 2	'F7'	274; ...
1	'A6'	290; 1	'D6'	290; 1	'B4'	323; 2	'B6'	323; 2	'A2'	355; ...
2	'E10'	484; 3	'B6'	581.0; 3	'F7'	612.5; 3	'B9'	645.0; 3	'B5'	645.0; ...
3	'A10'	371.0; 3	'F7'	452.0; 3	'F9'	475.5; 3	'B9'	306.0; 3	'B10'	217.0};
cmap = [.3 .3 .3; .5 .5 .5; .7 .7 .7];

for s = 1:3
idx = find([tdata{:, 1}]== s);
 counts(s,:) = histc(cat(1, tdata{idx, 3}), [0:100:700]);
end

figure(10); clf
b = bar([0:100:700],counts,'stacked');
for s = 1:3
    set(b(s), 'FaceColor', cmap(s, :));
end

val = []
for s = [1 3]
idx = find([tdata{:, 1}]== s);
 val = cat(1, val,cat(1, tdata{idx, 3}) );
end
round(prctile(val, [15 85]))






