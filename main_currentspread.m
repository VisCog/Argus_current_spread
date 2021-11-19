close  all
clear all

% electrode parameters
p.x1 = -2600; p.x2 = 2600; % separation of the electodes, max 
p.x1 = -1000; p.x2 = 1000;
p.lim = [2000, 900, 1100]; % how much retina to simulate
ss = 50; % simulation resolution

p.y = 0; 
p.rad = 225/2;

p.k = 6.5; % Ahuja params
p.k = [3.5 ]; % Esler
p.gamma = 1.69;
p.threshOrig = 23.392 *8.5; % 23 is spiking threshold, 50 a reasonable psycho estimate for a pulse train

% plot params
p.plot = 1;
p.pos = [ 42   1088  939   231];

% how bright the current is on the graphs
p.Isc = .4;

xv = -p.lim(1):ss:p.lim(1);
yv = -p.lim(2):ss:p.lim(2);
zv = -p.lim(3):ss:p.lim(3);
[p.X,p.Y,p.Z] = meshgrid(xv,yv,zv);

% find overlap function
ct = 1;
zvals = [500:100:3000]; % different heights
for z = zvals
    p.thresh = p.threshOrig; % each electrode goes up to 127
    p.z = -z;
    p.a = 400;
    
    p.R1 = sqrt(((p.X-p.x1).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
    p.R1 = p.R1-p.rad; p.R1(p.R1<0) = 0; p.R1 = p.R1./(max(p.R1(:)));
    p.R2 = sqrt(((p.X-p.x2).^2) + ((p.Y-p.y).^2)+((p.Z-p.z).^2));
    p.R2 = p.R2-p.rad; p.R2(p.R2<0) = 0; p.R2 = p.R2./(max(p.R2(:)));
    
    % find the threshold for the z
    [p,err] = fit('fit_currentspread',p, {'a'});
    out.amp(ct) = p.a;
    
    % now double it and show what the current spread looks like
    p.a = p.a*2;
    p.thresh = p.thresh*2;
    [I, p, Vq] = create_currentspread(p);   
    p.pos = [ 42   1088-((ct-1)*50)   939   231];
    figure(ct); clf
    create_currentspreadfig(I, p); drawnow;
    
    % current value directly under the electrode on retinal surface
    out.underelectrode(ct) = interpn(unique(p.Y), unique(p.X), unique(p.Z), I, 0, p.x2, 0);
    % max current value on the surface of the retina
    out.max(ct) = max(max(I(:, :, round(size(I, 3)/2))));
    % current value directly in between the two electrodes
    out.midpoint(ct) = interpn(unique(p.Y), unique(p.X), unique(p.Z), I, 0, 0, 0);
    set(gcf, 'Name', ['height ', num2str(z)]);
    ct = ct+1;
end

figure(100)
plot(zvals, out.amp, 'k')
xlabel('height from the retinal surface')
ylabel('threshold current')

figure(101); 
h(1) = plot(-zvals, out.underelectrode, 'r'); hold on
h(2) = plot(-zvals,out.midpoint, 'b'); hold on
h(3) = plot(-zvals,out.max, 'g--'); hold on
legend(h, {'under electrode', 'midpoint', 'max'});
xlabel('z'); ylabel('max current')


