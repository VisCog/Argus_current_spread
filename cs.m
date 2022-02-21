% cs.m
%
% holds all support functions for current spread
% project.
%
% functions can be called from outside with 'cs.<function name>'


classdef cs
    methods(Static)
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
                ret.e(i).I = ret.t./(1+((ret.k.*ret.e(i).R).^ret.a));
                ret.I = ret.I + ret.e(i).I;
            end
            ret.I(ret.Z<max([ret.z ret.z])) = 0; % remove all current above the electrode
            Sxy = ret.I(:, :, round(size(ret.I, 3)/2));
            if isfield(ret, 'loc') % if choosing a particular location on the retina at which to know amplitude
                disp('calculating threshold at a specific location')
                ret.Vq = interpn(unique(ret.Y), unique(ret.X),Sxy, ret.loc(2), ret.loc(1));
            else % find the max amplitude
                ret.Vq = max(Sxy(:));
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
                ret.t = (hi+lo)/2;
                ret = cs.create_currentspread(ret);
                if ret.t_ret < ret.Vq
                    hi = ret.t;
                else
                    lo = ret.t;
                end

                % if the fitting doesn't seem to work, check that the original hi and
                % lo ranges aren't too narrow
            end
            ret.t = (hi+lo)/2;
             disp(['max current on retina = ', num2str(ret.Vq), ' retinal thresh = ', num2str(ret.t_ret)]);
        end


        function I = create_currentspreadfig(ret)
            Sxy = ret.I(:, :, round(size(ret.I, 3)/2)); % slice on the retina, z = 0

            subplot(1,2,1)
            imagesc(unique(ret.X(:)), unique(ret.Y(:)), Sxy*ret.Isc); colormap(hot); hold on
            for i = 1:length(ret.x)
                h = viscircles([ret.x(i), ret.y(i)],ret.rad-1, 'Color', [.8 .8 .3], 'LineWidth', .5);
            end
            if isfield(ret, 'loc')
                x = unique(ret.X(:)); y = unique(ret.Y(:));
                for l = 1:size(ret.loc, 1)
                    plot(x(ret.loc(l, 1)),y(ret.loc(l, 2)),'g*');
                end
                plot(0, 0, 'b*')
            end
            if isfield(ret, 't_ret')
                contour(unique(ret.X(:)), unique(ret.Y(:)), Sxy, [-1 ret.t_ret*.95], 'g--');
            end
            axis equal
            title('retinal plane');

            subplot(1,2,2)
            Sxz = flipud(rot90(squeeze(ret.I(round(size(ret.I, 1)/2), :, :)))); % slice along the midpoint of y
            imagesc(unique(ret.X(:)), unique(ret.Z(:)),Sxz*ret.Isc); colormap(hot); hold on
            plot(unique(ret.X(:)), zeros(size(unique(ret.X(:)))), 'w--');
            for i = 1:length(ret.x)
                plot([ret.x(i)-ret.rad, ret.x(i)+ret.rad], [ret.z-3, ret.z-3],'-', 'Color', [.8 .8 .3],  'LineWidth', 2);
            end
            if isfield(ret, 't_ret')
                contour(unique(ret.X(:)), unique(ret.Z(:)), Sxz, [-1 ret.t_ret]*.95, 'g--');
            end
            title('depth plane');

            axis equal
            if isfield(ret, 'loc')
                for l = 1:size(ret.loc, 1)
                    plot(x(ret.loc(l, 1)), 0, 'g*');
                end
                  plot(0, 0, 'b*')
            end
        end
    end
end
