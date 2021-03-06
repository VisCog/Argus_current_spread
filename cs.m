% cs.m
%
% holds all support functions for current spread
% project.
%
% functions can be called from outside with 'cs.<function name>'


classdef cs
    methods(Static)
        function ret = setdefaultparams(flag)

            ret.rad = 225/2;
            ret.lim = [500, 1500, 1100]; % how much retina to simulate. x, y, z
            ret.Isc = .439; % how bright the current is on the graphs
            ret.eI = 400; % current on the electrode
            ret.t_ret_min = 50; % minimum retinal current required to reach perceptual threshold
            ret.t_ret = 50;
            if strcmp(flag.res, 'lowres')  % lower sampling for faster runtime
                ret.z_range = 0:200:2000;
                ret.d_range = 250:500:3000;
                ret.ss = 50; % simulation resolution
                ret.a_range = linspace(.5, 3, 5);
                ret.k_range = linspace(3, 15, 5); % controls current spread
                ret.rd_range = 1:5:25; % scaling factor representing retinal damage
            else
                % neurophysiological paramers
                ret.ss = 50;
                ret.z_range = 0:25:1300;
                ret.d_range = 250:50:3000;
                ret.a_range = .5:.25:3;
                ret.k_range = 3:.25:15; % controls current spread
                ret.rd_range = 1:.25:12; % threshold elevation due to retinal damage (mAmps)
            end
            xv = -ret.lim(1):ret.ss:ret.lim(1);
            yv = -ret.lim(2):ret.ss:ret.lim(2);
            zv = -ret.lim(3):ret.ss:ret.lim(3);
            [ret.X,ret.Y,ret.Z] = meshgrid(yv,xv,zv); % simulates the retina as a grid
        end

        function ret = calc_dist_from_electrode(ret)
            for i = 1:length(ret.x)
                ret.e(i).R = sqrt(((ret.X-ret.x(i)).^2) + ((ret.Y-ret.y(i)).^2)+((ret.Z-ret.z).^2));
                ret.e(i).R = ret.e(i).R-ret.rad;
                ret.e(i).R(ret.e(i).R<0) = 0;
                ret.e(i).R = ret.e(i).R./(max(ret.e(i).R(:))); % deal with directly under the electrode
            end
        end
        function ret = create_currentspread(ret)
            % takes in parameter values describing electrodes and
            % calculates the resulting current field

            ret.I = zeros(size(ret.e(1).R));
            for i = 1:length(ret.x) % for each electrode
                ret.e(i).I = ret.eI./(1+((ret.k.*ret.e(i).R).^ret.a));
                ret.I = ret.I + ret.e(i).I;
            end
            ret.I(ret.Z<max([ret.z ret.z])) = 0; % remove all current above the electrode
            eIxy = ret.I(:, :, round(size(ret.I, 3)/2));
            if isfield(ret, 'loc') & ~isempty(ret.loc) % if choosing a particular location on the retina at which to know amplitude
                disp('calculating threshold at a specific location')
                ret.eIxy = interpn(unique(ret.Y), unique(ret.X), eIxy, ret.loc(2), ret.loc(1));
            else % find the max amplitude
                ret.eIxy = max(eIxy(:));
            end
        end

        function [err] = fit_currentspread(ret)
            % calculates current spread, and finds the electrode amplitude
            % which produces a current of ret.t_ret on the retinal surface
            ret = cs.create_currentspread(ret);
            err = (ret.Vq(1)-ret.t_ret).^2;
        end

        function ret = fit_currentspreadfast(ret, fitParams)

            lo = fitParams.lo; hi = fitParams.hi;
            for i=1:fitParams.nreps
                ret.eI = (hi+lo)/2;
                ret = cs.create_currentspread(ret);
                if ret.t_ret < ret.eIxy
                    hi = ret.eI;
                else
                    lo = ret.eI;
                end

                % if the fitting doesn't seem to work, check that the original hi and
                % lo ranges aren't too narrow
            end
            ret.eI = (hi+lo)/2;
            % disp(['max current on retina = ', num2str(ret.eIxy), ' retinal thresh = ', num2str(ret.t_ret)]);
        end
        function ret = calculate_dip(ret, varargin)
            % second optional input:
            % flag.dip =    'perceptual', default is assuming that electric current below threshold has no perceptual effect
            %               'electrical', calculates dip based on electric
            % field, no effect of retinal threshold
            % flag.rd does the effect of rd just change the 'threshold' or
            % 'scale' the effectiveness of current by a constant factor
    
            if nargin <2
                flag.dip = 'perceptual'; % default is assuming that electric current below threshold has no perceptual effect, alternative is 'electrical'
                flag.rd = 'scale';
            else
                flag = varargin{1};
            end
            if strcmp(flag.rd, 'scale') % the scaling takes account of the change in threshold
                Itmp = ret.I./2;
                rmin = ret.t_ret_min;
            else
                Itmp = ret.I;
                rmin = ret.t_ret; % include rd int the threshold
            end
            Ixy = Itmp(:, :, round(size(Itmp, 3)/2));
            % max current value on the surface of the retina
            ret.I_max = max(Ixy(:));
            ret.I_mid = interpn(unique(ret.Y), unique(ret.X), unique(ret.Z), Itmp, 0, 0, 0); % current value directly in between the two electrodes


            ret.I_max = ret.I_max - rmin; ret.I_max(ret.I_max<0)=0;
            ret.I_mid = ret.I_mid - rmin; ret.I_mid(ret.I_mid<0)=0;

            ret.dip = 100.*(ret.I_max-ret.I_mid)./ret.I_max;
            % location of the max
            [y, x] = find(Ixy == ret.I_max);
            ret.loc = [x y];
        end


        function I = create_currentspreadfig(ret)
            Sxy = ret.I(:, :, round(size(ret.I, 3)/2)); % slice on the retina, z = 0
            disp(['max I',  num2str(round(max(ret.I(:))))]);
            subplot(1,2,1)
            image(unique(ret.X(:)), unique(ret.Y(:)), Sxy*ret.Isc); colormap(hot(256)); hold on
            for i = 1:length(ret.x)
                h = viscircles([ret.x(i), ret.y(i)],ret.rad-1, 'Color', [.8 .8 .3], 'LineWidth', .5);
            end
            if isfield(ret, 'loc')
                x = unique(ret.X(:)); y = unique(ret.Y(:));
                for l = 1:size(ret.loc, 1)
                    plot(x(ret.loc(l, 1)),y(ret.loc(l, 2)),'bp', 'MarkerSize', 5);
                end
                plot(0, 0, 'cp', 'MarkerSize', 5)
            end
            if isfield(ret, 't_ret')
                contour(unique(ret.X(:)), unique(ret.Y(:)), Sxy, [-1 ret.t_ret*.95], 'b--', 'LineWidth', 2);
            end
            axis equal
            title('retinal plane');

            subplot(1,2,2)
            Sxz = flipud(rot90(squeeze(ret.I(round(size(ret.I, 1)/2), :, :)))); % slice along the midpoint of y
            disp(['max I ',  num2str(round(max(ret.I(:))))]);
            image(unique(ret.X(:)), unique(ret.Z(:)),Sxz*ret.Isc); 
            colormap(hot(256)); hold on
           % colorbar(gca)
            plot(unique(ret.X(:)), zeros(size(unique(ret.X(:)))), 'w--');
            for i = 1:length(ret.x)
                plot([ret.x(i)-ret.rad, ret.x(i)+ret.rad], [ret.z-3, ret.z-3],'-', 'Color', [.8 .8 .3],  'LineWidth', 2);
            end
            if isfield(ret, 't_ret')
                contour(unique(ret.X(:)), unique(ret.Z(:)), Sxz, [-1 ret.t_ret]*.95, 'b--', 'LineWidth', 2);
            end
            title('depth plane');

            axis equal
            if isfield(ret, 'loc')
                for l = 1:size(ret.loc, 1)
                    plot(x(ret.loc(l, 1)), 0, 'bp', 'MarkerSize', 5);
                end
                plot(0, 0, 'b*')
            end        
        end

        function [xc,yc] = unwrap_contour(c, cval)
            count = 1;
            xc = [];
            yc = [];
            i=1;

            while count<size(c,2)
                cc = c(1,count);% contour level
                n = c(2,count); % number of vertices in that contour level
                if cc == cval
                    xc = c(1,(count+1):(count+n));
                    yc = c(2,(count+1):(count+n));
                end

                count = count+n+1;
                i=i+1;
            end
        end
        function prctile_out = plot_amp_thresholds(prctile_vals)
            % plot thresholds and return amplitudes within a certain pctile
            %% thresholds
            tdata = {1	'B10'	89.0; 1	'B10'	93.0; 2	'A8'	153; 1	'B9'	177; ...
                1	'F10'	177; 2	'A4'	194; 1	'F10'	202; 1	'C6'	210; 2	'F2'	226; ...
                1	'E9'	242; 1	'A8'	250; 2	'D1'	274; 2	'F2'	274; 2	'F7'	274; ...
                1	'A6'	290; 1	'D6'	290; 1	'B4'	323; 2	'B6'	323; 2	'A2'	355; ...
                2	'E10'	484; 3	'B6'	581.0; 3	'F7'	612.5; 3	'B9'	645.0; 3	'B5'	645.0; ...
                3	'A10'	371.0; 3	'F7'	452.0; 3	'F9'	475.5; 3	'B9'	306.0; 3	'B10'	217.0};
            cmap = [ 0 .6 .4; ; 0 0 1; 0 .5 1 ];

            % calculate and plot amplitude values
            for s = 1:3
                idx = find([tdata{:, 1}]== s);
                counts(s,:) = histc(cat(1, tdata{idx, 3}), [0:100:700]);
            end

            figure(10); clf
           b = bar([0:100:700],counts,'stacked', 'EdgeAlpha', 0);

            for s = 1:3
                set(b(s), 'FaceColor', cmap(s, :));
            end

            % collate values across subjects
            val = [];
            for s = [1 2 3]
                idx = find([tdata{:, 1}]== s);
                val = cat(1, val,cat(1, tdata{idx, 3}) );
            end

            disp(['reporting percentiles across all subjects' num2str(prctile_vals)]);
            prctile_out = round(prctile(val, prctile_vals));
            disp(prctile_out);

            disp('reporting mean and percentiles for each subject')
            val = [];
            for s = [1 2 3]
                idx = find([tdata{:, 1}]== s);
                disp(['subject', num2str(s)])
                disp(['number of electrodes  = ', num2str(length(idx))]);
                disp(['mean  = ', num2str(round(mean(cat(1, tdata{idx, 3})), 1))]);
                disp(['percentiles  = ', num2str(round(prctile(cat(1, tdata{idx, 3}), prctile_vals),1))]);
            end
        end
    end
end

