function I = create_currentspreadfig(I,p)

I(p.Z<max([p.z p.z])) = 0;

Sxy = I(:, :, round(size(I, 3)/2)); % slice on the retina, z = 0
subplot(1,2,1)
image(unique(p.X(:)), unique(p.Y(:)), Sxy*p.Isc); colormap(hot); hold on
h = viscircles([p.x1, p.y],p.rad-1, 'Color', [.8 .8 .3], 'LineWidth', .5);
h = viscircles([p.x2, p.y],p.rad-1, 'Color', [.8 .8 .3], 'LineWidth', .5);
if isfield(p, 'loc')
plot(p.loc(1), p.loc(2), 'g*');
end
contour(unique(p.X(:)), unique(p.Y(:)), Sxy, [-1 p.thresh*.95], 'g--');
axis equal

title('retinal plane');
% axis off
set(gcf, 'Position', p.pos);
subplot(1,2,2)
Sxz = flipud(rot90(squeeze(I(round(size(I, 1)/2), :, :)))); % slice along the midpoint of y
image(unique(p.X(:)), unique(p.Z(:)),Sxz*p.Isc); colormap(hot); hold on
plot(unique(p.X(:)), zeros(size(unique(p.X(:)))), 'w--');
plot([p.x1-p.rad, p.x1+p.rad], [p.z-3, p.z-3],'-', 'Color', [.8 .8 .3],  'LineWidth', 2);
plot([p.x2-p.rad, p.x2+p.rad], [p.z-3, p.z-3], '-', 'Color', [.8 .8 .3],   'LineWidth', 2);
contour(unique(p.X(:)), unique(p.Z(:)), Sxz, [-1 p.thresh*0.95], 'g--');
title('depth plane');
set(gcf, 'Position',p.pos );

axis equal
if isfield(p, 'loc')
plot(p.loc(1), p.loc(3), 'g*');
end
end
