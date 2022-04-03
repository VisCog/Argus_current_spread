% code for Yucel et al
%
% to replicate results for initial submission select:
% flag.res = highres;
% flag.dip = perceptual;
% flag.rd = 'threshold'
%
% but better parameters are probably
% flag.res = highres;
% flag.dip = perceptual;
% flag.rd = 'scale'

clear all

%% plot amplitude thresholds
figure(1); clf
prctile_out = cs.plot_amp_thresholds([15 50 85]); % get percentile range for amplitude thresholds
% fitting params
fitParams.nreps = 12; % shouldn't be below 10
fitParams.tol = 0.1; fitParams.lo = 0; fitParams.hi = 5000; fitParams.thr = 1;

%% begin simulations, basic parameters
safety_lim = 660;
flag.res = 'highres'; % running the code with coarse or fine parameter sampling. Warning! highres is very slow!
flag.dip = 'perceptual';  % default is assuming that electric current below threshold has no perceptual effect, alternative is 'electrical'
flag.rd = 'scale'; % 'threshold', RD just changes the threshold, or scales all current linearly
savestr = ['4_2_2022', '_', flag.res];


%% calculate thresholds as a function of height, for all the possible retinal damages, single electrode
ret = cs.setdefaultparams(flag);
ret.x = 0; ret.y = 0; % location of the electrode
ct = 1;
for z = 1:length(ret.z_range)
    disp(['fitting z = ', num2str(z),  ' out of ', num2str(length(ret.z_range))]);
    ret.z = -ret.z_range(z);
    ret = cs.calc_dist_from_electrode(ret);
    for r = 1:length(ret.rd_range)
        ret.t_ret = ret.rd_range(r) * ret.t_ret_min;
        for a = 1:length(ret.a_range)
            ret.a = ret.a_range(a);
            for ks = 1:length(ret.k_range)
                ret.k = ret.k_range(ks);
                % find the threshold for the z
                ret = cs.fit_currentspreadfast(ret, fitParams);
                if ret.eI <= safety_lim % if this threshold below the safety limit
                    ret_single.eI(ct) = ret.eI;
                    ret_single.a(ct) = ret.a_range(a);   ret_single.k(ct) =ret.k_range(ks);
                    ret_single.z(ct) = ret.z_range(z);
                    ret_single.Th_RD(ct) =ret.rd_range(r); % how much retinal damage raises threshold
                    ret_single.t_ret(ct) = ret.t_ret;
                    if z > 1
                        ind = find(ret_single.a ==ret.a_range(a) & ret_single.k==ret.k_range(ks) & ...
                            ret_single.Th_RD==ret.rd_range(1) &  ret_single.z==ret.z_range(1));
                        ret_single.Th_Z(ct) = ret_single.eI(ct)./(ret_single.eI(ind) * ret_single.Th_RD(ct)); % calculate how much lift is raising threshold
                    end
                    ct = ct + 1;
                end; end; end; end; end
disp([num2str(ct-1), ' parameterizations have thresholds below the safety limit']);


%% plot thresholds as a function of height, for a variety of retinal damage, single electrode
figure(1); clf

for a = 1:length(ret.a_range)
    for ks = 1:length(ret.k_range)
        for r = 1:length(ret.rd_range)
            ind = find(ret_single.a ==ret.a_range(a) & ret_single.k==ret.k_range(ks) & ...
                ret_single.Th_RD==ret.rd_range(1));

            subplot(1,2,1)
            plot(ret_single.z(ind), ret_single.eI(ind), 'k'); hold on
            xlabel('z'); ylabel('threshold')
            set(gca, 'XLim', [0 1000]); set(gca, 'XTick', [0:200:1000])
            set(gca, 'YLim', [0 700]);set(gca, 'YTick',[0:100:700])

            subplot(1,2,2)
            plot(log(ret_single.z(ind)), log(ret_single.eI(ind)) , 'k'); hold on
            xlabel('z'); ylabel('threshold non damaged retina')
            set(gca, 'XTick', log([ 100 1000]));  set(gca, 'XLim', log([50 2000])); set(gca, 'XTickLabel',[100 1000]); logx2raw
            set(gca, 'YTick', log([10  100 1000])); set(gca, 'YLim', log([9 2000])); set(gca, 'YTickLabel',[10 100 1000]); logy2raw
        end
    end
end
save([savestr, 'single']);

%% plot, current spreads for 2 electrodes at different heights
fitParams.nreps = 20;
clear ret; ret = cs.setdefaultparams(flag);
ret.ss = 25;
ret.a = 1.5; ret.k = 15;
ret.t_lift = 0;
ret.x = [-700 700];
ret.y = [0 0];
ret_saved = ret;
plt.pos = [ 42   1088  939   231];
z_list = [-150 -750];
for zz = 1:2 % simulating pairs of electrodes, at two different heights
    ret.z = z_list(zz);
    figure(zz+1);clf
    ret = cs.calc_dist_from_electrode(ret);
    if isfield(ret, "loc")
        ret = rmfield(ret, "loc");
    end
    ret = cs.fit_currentspreadfast(ret, fitParams);
    % now double it and calculate what the current spread looks like
    ret.eI =  ret.eI * 2;
    ret = cs.create_currentspread(ret);
    ret = cs.calculate_dip(ret, flag);

    cs.create_currentspreadfig(ret); set(gcf, 'Position', plt.pos);
    set(gcf, 'Name', ['height = ', num2str(ret.z)]); drawnow
end


%% simulate a range of parameter values
clear ret; ret = cs.setdefaultparams(flag);
ret.y = [0 0];
ct = 1;
for d = 1:length(ret.d_range)
    disp(['fitting d = ', num2str(d),  ' out of ', num2str(length(ret.d_range))]);
    ret.x = [-round(ret.d_range(d)/2) round(ret.d_range(d)/2)];
    for s = 1:length(ret_single.eI) % for each parameterization below safety limits
        ret.z = -ret_single.z(s);
        ret.t_ret = ret_single.t_ret(s);
        ret.k = ret_single.k(s);
        ret.a = ret_single.a(s);
        ret.eI = ret_single.eI(s) * 2;
        if ret.eI>safety_lim
            ret.eI = safety_lim;
        end
        ret = cs.calc_dist_from_electrode(ret);
        ret = cs.create_currentspread(ret);
        ret = cs.calculate_dip(ret, flag);

        ret_pair.eI(ct) = ret_single.eI(s);
        ret_pair.a(ct) = ret_single.a(s);
        ret_pair.k(ct) = ret_single.k(s);
        ret_pair.z(ct) = ret_single.z(s);
        ret_pair.Th_RD(ct) = ret_single.Th_RD(s);
        ret_pair.Th_Z(ct) = ret_single.Th_Z(s);
        ret_pair.ret.t_ret(ct) = ret_single.t_ret(s);
        ret_pair.dist(ct) = ret.d_range(d);
        ret_pair.I_max(ct)= ret.I_max;
        ret_pair.I_mid(ct) = ret.I_mid;
        ret_pair.dip(ct) = ret.dip;

        ct = ct+1;
    end
end
save([savestr, '_pair']);

%% create isodipcurve
% replot data from the regression
figure(4); clf; hold on


p = 1.0e+02 * [0.003063282575664  -1.546555906566781];
amp_pred = polyval(p, ret.d_range);
plot(ret.d_range(2:end),amp_pred(2:end) , 'g-', 'LineWidth', 2); hold on
set(gca, 'YLim', [0 700])
set(gca, 'XLim', [500 4000])
xlabel('Physical Distance ');
ylabel('Amplitude');

amp_max = 499; %280; 177; %;
amp_min = 177; %177;%177; %

% find the distances corresponding to each possible dip criterion, and the
% expected amplitude based on the regression model, for that distance
clear keepSim;
err_Thr = 50;
ct = 1;
for ks = 1:length(ret.k_range)
    for a = 1:length(ret.a_range)
        for r = 1:length(ret.rd_range)
            clear M;
            ind = find(ret_pair.eI>amp_min & ret_pair.eI<amp_max & ret_pair.a == ret.a_range(a) & ret_pair.k == ret.k_range(ks)  & ret_pair.Th_RD == ret.rd_range(r));
            dip_min = max([ceil(min(ret_pair.dip(ind))), 10]);
            dip_max = min([floor(max(ret_pair.dip(ind))), 90]);
            % for a given dip'
            Th_Z_unique = unique(ret_pair.Th_Z(ind));
            dist_unique = unique(ret_pair.dist(ind));
            dip_vals = ret_pair.dip(ind);
            Th_Z_vals = ret_pair.Th_Z(ind);
            z_vals = ret_pair.z(ind);
            Th_RD_vals = ret_pair.Th_RD(ind);
            eI_vals = ret_pair.eI(ind);
            dist_vals = ret_pair.dist(ind);

            % calculate the matrix of dip as a function of inter-electrode
            % distance and Z. This matrix doesn't vary with RD
            if length(Th_Z_unique)>2 & length(dist_unique)>2
                for i = 1:length(Th_Z_unique)
                    for j = 1:length(dist_unique)
                        gval = find(Th_Z_vals == Th_Z_unique(i) & dist_vals == dist_unique(j) & Th_RD_vals == ret.rd_range(r));
                        M(j, i) = nanmean(dip_vals(gval));
                    end;  end

                [c,h] =contour(Th_Z_unique,dist_unique, M, [dip_min:dip_max]);
                for dd = dip_min:dip_max

                    [Th_Z_c, dist_c] = cs.unwrap_contour(c, dd); % the Z values as a function of separation for a single iso curve
                    eI_c = interp1(Th_Z_vals(1:length(Th_Z_unique)), eI_vals(1:length(Th_Z_unique)), Th_Z_c);
                    z_c = interp1(Th_Z_vals(1:length(Th_Z_unique)), z_vals(1:length(Th_Z_unique)), Th_Z_c);

                    g_eI = find(eI_c>amp_min & eI_c<amp_max);

                    pred_eI_c = polyval(p, dist_c(g_eI));
                    err = sqrt(mean((eI_c(g_eI)-pred_eI_c).^2));
                    if err<err_Thr
                        keepSim(ct).eI = eI_c(g_eI);
                        keepSim(ct).dist = dist_c(g_eI);
                        keepSim(ct).Th_Z = Th_Z_c(g_eI);
                        keepSim(ct).z = z_c(g_eI);
                        keepSim(ct).Th_RD = ones(size(Th_Z_c(g_eI))).*ret.rd_range(r);

                        keepSim(ct).err = err;
                        keepSim(ct).a = ret.a_range(a);
                        keepSim(ct).k = ret.k_range(ks);
                        keepSim(ct).dip = dd;
                        [~, si] = sort(keepSim(ct).dist); % probably unnecessary
                        patchline( keepSim(ct).dist(si), keepSim(ct).eI(si), ...
                            'EdgeColor',[.3 .3 .3], 'LineWidth',2,'EdgeAlpha',0.2); drawnow
                        ct = ct+1;
                    end; end;  end;  end;  end; end


%% now calculate what would have happened without lift, RD or either
clear ret; ret = cs.setdefaultparams(flag);

for ks = 1:length(keepSim)
    ret.a = keepSim(ks).a;
    ret.k = keepSim(ks).k;
    ret.y = [0 0];
    for ps = 1:length(keepSim(ks).eI)
        ret.eI =  keepSim(ks).eI(ps) * 2;
        if ret.eI>safety_lim
            ret.eI = safety_lim;
        end

        % no RD
        ret.z = -keepSim(ks).z(ps);
        ret.t_ret  = ret.t_ret_min;
        lo = 1; hi = 3000;
        for i = 1:10
            dist = (hi+lo)/2;
            ret.x = [-round(dist/2) round(dist/2)];
            ret = cs.calc_dist_from_electrode(ret);
            ret = cs.create_currentspread(ret);
            ret = cs.calculate_dip(ret, flag);

            if ret.dip >  keepSim(ks).dip
                hi = dist;
            else
                lo = dist;
            end
        end
        keepSim(ks).no_RD_dist(ps) = (hi+lo)/2;
        keepSim(ks).no_RD_dip(ps) = ret.dip;

        % no Z
        ret.z = 0;
        ret.t_ret  = ret.t_ret_min * keepSim(ks).Th_RD(ps);
        lo = 1; hi = 3000;
        for i = 1:10
            dist = (hi+lo)/2;
            ret.x = [-round(dist/2) round(dist/2)];
            ret = cs.calc_dist_from_electrode(ret);
            ret = cs.create_currentspread(ret);
            ret = cs.calculate_dip(ret, flag);
            if ret.dip >  keepSim(ks).dip
                hi = dist;
            else
                lo = dist;
            end
        end
        keepSim(ks).no_Z_dist(ps) = (hi+lo)/2;
        keepSim(ks).no_Z_dip(ps) = ret.dip;

        % no RD or Z
        ret.z =  0 ;
        ret.t_ret  = ret.t_ret_min;
        lo = 1; hi = 3000;
        for i = 1:10
            dist = (hi+lo)/2;
            ret.x = [-round(dist/2) round(dist/2)];
            ret = cs.calc_dist_from_electrode(ret);
            ret = cs.create_currentspread(ret);
            ret =  cs.calculate_dip(ret, flag);
            if ret.dip >  keepSim(ks).dip
                hi = dist;
            else
                lo = dist;
            end
        end
        keepSim(ks).no_Z_RD_dist(ps) = (hi+lo)/2;
        keepSim(ks).no_Z_RD_dip(ps) = ret.dip;
    end
end

%% scatter plots

%% Z vs. RD for successful simulations
figure(5); clf
for ks = 1:length(keepSim)
    p = scatter([keepSim(ks).Th_Z], [keepSim(ks).Th_RD]+.03*randn(size(keepSim(ks).Th_RD)), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .5); hold on
end
xlabel('Lift values');
ylabel('RD values');
axis equal
set(gca, 'XLim', [.9 9])
set(gca, 'YLim', [.9 5.5])

%% Histogram of dip values for successful simulations
figure(6); clf
hist([keepSim(:).dip])
set(gca, 'XLim', [0 100])
set(gca, 'XTick', 0:10:100)
xlabel('Dip required for 60% discrimination')
set(gcf, 'Position', [1000        1107        1121         231])

%% histogram of distances, when Z, RD or both is removed
figure(7); clf
orig = []; no_Z = []; no_RD = []; no_Z_RD = [];
for ks = 1:length(keepSim)
    orig = cat(1, orig, keepSim(ks).dist(:));
    no_Z = cat(1, no_Z, keepSim(ks).no_Z_dist(:));
    no_RD = cat(1, no_RD, keepSim(ks).no_RD_dist(:));
    no_Z_RD = cat(1, no_Z_RD, keepSim(ks).no_Z_RD_dist(:));
end
ss = 50;
h = histogram(orig,100:ss:3000, 'FaceColor',[1 1 1], 'FaceAlpha', .3, 'EdgeAlpha', 1); hold on
h_no_Z = histogram(no_Z,100:ss:3000, 'FaceColor',[1 0 0], 'FaceAlpha', .3, 'EdgeAlpha',1); hold on
h_no_RD = histogram(no_RD,100:ss:3000, 'FaceColor',[0 1 0], 'FaceAlpha', .3, 'EdgeAlpha',1); hold on
h_no_Z_RD = histogram(no_Z_RD,100:ss:3000, 'FaceColor',[0 0  1], 'FaceAlpha', .3, 'EdgeAlpha', 1); hold on

figure(8);clf
plot(h.BinEdges(2:end)-h.BinWidth/2, h.Values./sum(h.Values), 'k', 'LineWidth', 2); hold on
plot(h_no_Z.BinEdges(2:end)-h_no_Z.BinWidth/2, h_no_Z.Values./sum(h_no_Z.Values), 'r','LineWidth', 2); hold on
plot(h_no_Z_RD.BinEdges(2:end)-h_no_Z_RD.BinWidth/2, h_no_Z_RD.Values./sum(h_no_Z_RD.Values),  'Color', [.5 0 1], 'LineWidth', 2); hold on
plot(h_no_RD.BinEdges(2:end)-h_no_RD.BinWidth/2, h_no_RD.Values./sum(h_no_RD.Values), 'b','LineWidth', 2); hold on

line([1557 1557], [ 0 .35],'Color', [.5 .5 .5]); hold on
line([2291 2291], [ 0 .35], 'Color', [.5 .5 .5])
line([1324 1324], [ 0 .35], 'Color', [.5 .5 .5])
%set(gca, 'YLim', [ 0 45])
line([1557 1557], [ 0 .35],'Color', [.5 .5 .5]); hold on
line([2291 2291], [ 0 .35], 'Color', [.5 .5 .5])
line([1324 1324], [ 0 .35], 'Color', [.5 .5 .5])
set(gca, 'XLim', [0 2700])
set(gca, 'YLim', [0 .35])

xlabel('Min Distance required for 60% discrimination, pctiles')
disp(['distances all :', num2str(round(prctile(orig, [25 50 70])))]); hold on
disp(['no z :', num2str(round(prctile(no_Z, [25 50 70])))]); hold on
disp(['no rd :', num2str(round(prctile(no_RD, [25 50 70])))]); hold on
disp(['no z or rd :', num2str(round(prctile(no_Z_RD, [25 50 70])))]); hold on


