close  all


% electrode parameters
p.x1 = 0; p.x2 = 0; % separation of the electodes, max 
p.lim = [1000, 200, 400]; % how much retina to simulate
ss = 10; % simulation resolution

p.y = 0; p.z = -100;
p.rad = 15;
p.a = 12.09;

p.k = [3.5 ]; % esler
%p.k = 6.5;  % Ahuja params
p.gamma = 1.69;

% plot params
p.plot = 1;
p.pos = [ 42  1088  939   231];

% how bright the current is on the graphs
p.Isc = 10;

xv = -p.lim(1):ss:p.lim(1);
yv = -p.lim(2):ss:p.lim(2);
zv = -p.lim(3):ss:p.lim(3);

[p.X,p.Y,p.Z] = meshgrid(xv,yv,zv);

p.R1 = sqrt(((p.X-p.x1).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
p.R1 = p.R1-p.rad; p.R1(p.R1<0) = 0; p.R1 = p.R1./(max(p.R1(:)));
p.R2 = sqrt(((p.X-p.x2).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
p.R2 = p.R2-p.rad; p.R2(p.R2<0) = 0; p.R2 = p.R2./(max(p.R2(:)));



[I, p, Vq] = create_currentspread(p); 
p.pos = [ 42   1088   939   231];
figure(1); clf
create_currentspreadfig(I, p); drawnow;
colormap(flipud(hot));
%colormap((hot));


