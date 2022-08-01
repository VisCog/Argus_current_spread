
clear all
close all
fileName = 'all_subjects_summarized.csv';

% Load in the csv file
S = xls2struct(fileName);

%% List of free parameters for all models
funName = 'logit';  % doing logistic regression
% funName = 'normcdf';p

freeList_b0 = {'b0'};
freeList_a = {'b0','ka'};
freeList_dte = {'b0','kdte'};
freeList_dta = {'b0','kdta'};
freeList_a_dte = {'b0','ka','kdte'};
freeList_a_dta = {'b0','ka','kdta'};
freeList_dte_dta = {'b0','kdte','kdta'};
freeList_a_dte_dta = {'b0','ka','kdte','kdta'};

%% pull out the data into sensible variable names

sub = {S.subject_id{:}};  % subject  name
amp1 = [S.pts_amp1{:}];
amp2 = [S.pts_amp2{:}];
e1 ={S.pts_electrode1{:}}; % not used
e2 = {S.pts_electrode2{:}}; % not used
dte =  [S.electrode_distance{:}];  % distance to electrode
dta = [S.dta_bi{:}]; % distance to electrode - distance to the axon
daa = [S.daa_bi{:}];  % distance along the axon (not used)
a = (amp1+amp2)/2;  % mean amplitude

np = [S.pts_number_processed{:}];  % number of percepts reported (1 or 2)
prob2 = [S.prob_2{:}];

% list of subjects (not used since we're collapsing across subjects)
subs = unique(sub);
n = length(subs);  % number of subjects

%% Fitting and plotting for just amplitude and dte

% The probability of seeing 1 percept is defined by:
% 1. Generating a single 'y' axis as a linear combination of amplitude and
% distances, with constants p.ka and p.kdte.  This generates a plane in
% the amplitude-distance plane.
%
% 2. Passing this single value through a logistic function 1/(1+exp(-y))

% resp is the probability of 2 percepts
resp = prob2; % 1 percept: resp = 0, 2 percepts: resp =1

% initial parameters
% p.sd = 1;  % standard deviation of cumulative normal
% (kept constant at this large number so slopes are reasonable)

p.b0 = .5;
p.ka = 0;
p.kdte = 0;
p.kdta = 0;
p = fit('getErr',p,freeList_a_dte,dte,a,dta,resp,funName);
[err,x] =getErr(p,dte,a,dta,resp,funName);
[aic,bic] = aicbic(-err,2);

disp(['aic current spread = ', num2str(aic)]);
% Plot in 3D: psychometric function surface, and binned responses

% for the surface:
surfRez = 21;
binRez = 5;
%[dteList,aList] = meshgrid(linspace(min(dte),max(dte),surfRez),linspace(min(a),max(a),surfRez));
[dteList,aList] = meshgrid(linspace(575,max(dte),surfRez),linspace(50,max(a),surfRez));


dtaList = zeros(size(dteList));
%     xList =   p.ka*aList + p.kd*dList;
%     probList = normcdf(xList,0,p.sd);
[logL,~,probList] = getErr(p,dteList,aList,dtaList,[],funName);

ptmp = fit('getErr',p,freeList_a_dte_dta,dte,a,dta,resp,funName);
[err,x] =getErr(p,dte,a,dta,resp,funName);

[aic,bic] = aicbic(-err,3);
disp(['aic axon = ', num2str(aic)]);

figure(1)
clf
surf(dteList,aList,probList,'FaceAlpha',.5 ,'EdgeColor','none');
grid off
xlabel('Physical Distance (\mum) ');
ylabel('Amplitude (\muA)');
zlabel('P("2 Percepts")');
set(gca, 'XLim', [0 6000])
set(gca, 'YLim', [0 650])
% Plot data for each of the bins:
hold on

aBin = linspace(min(a),max(a),binRez);
dBin = linspace(min(dte),max(dte),binRez);

for i=1:(length(aBin)-1)
    for j=1:(length(dBin)-1)
        id = a>=aBin(i) & a<=aBin(i+1) & dte>=dBin(j) & dte<=dBin(j+1);
        probBin = mean(resp(id));
        if sum(id)>0
            plot3((dBin(j)+dBin(j+1))/2,(aBin(i)+aBin(i+1))/2,probBin,'ko',...
                'MarkerSize', sum(id)/4+10,'MarkerFaceColor',[1,1,1]*probBin);
        end
    end
end
c=  colorbar;
pos  = c.Position;
set(c,'Position',[0.820357142857141 0.19952380952381 0.031499999999999 0.493333333333336])
view(-20,20)
title(sprintf('b0 = %5.5f, ka = %5.5f, kd = %5.5f',p.b0,p.ka,p.kdte));

exportgraphics(gcf,['.' filesep 'figures' filesep 'fig3a.eps'],'ContentType','vector')
exportgraphics(gcf,['.' filesep 'figures' filesep 'fig3a.emf'],'ContentType','vector')
exportgraphics(gcf,['.' filesep 'figures' filesep 'fig3a.pdf'],'ContentType','vector')

%% Contour plot for amp and dte


contourList = [.65,.75,.85]; %two-point discrimination thresholds
f2 = figure(2);
f2.Position = [488 308 385 453];
clf
[c,h]= contour(dteList,aList,probList,contourList);

count = 1;
for i=1:length(contourList)
    n = c(2,count);
    xc{i} = c(1,(count+1):(count+n));
    yc{i} = c(2,(count+1):(count+n));
    count = count+n+1;
end

lineStyles = {'k:','k-','k-.'};
clf
hold on
for i=1:length(contourList)
    plot(xc{i},yc{i},lineStyles{i},'LineWidth',2)
    str{i} = sprintf('%d%%',100*contourList(i));
    text(xc{i}(end)-150,yc{i}(end)+25,str{i},'HorizontalAlignment','center');
    tmp = interp1(yc{i}, xc{i}, [210 274 476]);
    disp(['predicted spatial res = ', num2str(tmp)])
end
xlabel('Physical Distance (\mum)');
ylabel('Amplitude (\muA)');
text(1000,600,'One Percept');
text(2800,100,'Two Percepts');
%line([500 4000], [ 274 274], 'lineWidth',1.5)
%line([500 4000], [ 476 476], 'lineWidth',1.5)
%line([500 4000], [ 210 210], 'lineWidth',1.5)

%find 75% 2p thresh from median thresholds
tp_val = interp1(yc{2}, xc{2}, [210 274 476]);
disp(['predicted spatial res = ', num2str(tp_val)])
rectangle('Position',[0 0  tp_val(1) 210],'LineStyle','--', ...
          'EdgeColor', [0.5 0.5 0.5 0.5], 'LineWidth',2)
rectangle('Position',[0 0  tp_val(2) 274],'LineStyle','--', ...
          'EdgeColor', [0.5 0.5 0.5 0.5], 'LineWidth',2)
rectangle('Position',[0 0  tp_val(3) 476],'LineStyle','--', ...
          'EdgeColor', [0.5 0.5 0.5 0.5], 'LineWidth',2)
%legend(str,'Location','NorthWest')
set(gca,'XLim',1000*[.575,4])
set(gca,'YLim',[50,700])


exportgraphics(gcf,['.' filesep 'figures' filesep 'fig3b.eps'],'ContentType','vector')
exportgraphics(gcf,['.' filesep 'figures' filesep 'fig3b.emf'],'ContentType','vector')
exportgraphics(gcf,['.' filesep 'figures' filesep 'fig3b.pdf'],'ContentType','vector')
%% Fits, Chi-squared and p-values

% initial parameters
p.b0 = 0.5;  % constant
p.kdte = 0;  % slope constant for distance
p.kdta = 0;
p.ka =0;

% do all the model fits
p_b0 = fit('getErr',p,freeList_b0,dte,a,dta,resp,funName);
p_a = fit('getErr',p,freeList_a,dte,a,dta,resp,funName);
p_dte = fit('getErr',p,freeList_dte,dte,a,dta,resp,funName);
p_a_dte = fit('getErr',p,freeList_a_dte,dte,a,dta,resp,funName);
p_a_dta = fit('getErr',p,freeList_a_dta,dte,a,dta,resp,funName);
p_dte_dta = fit('getErr',p,freeList_dte_dta,dte,a,dta,resp,funName);
p_a_dte_dta = fit('getErr',p,freeList_a_dte_dta,dte,a,dta,resp,funName);

% get all the 'err's.  Err is the -log likelihood
[err_b0,x] =getErr(p_b0,dte,a,dta,resp,funName);
[err_a,x] =getErr(p_a,dte,a,dta,resp,funName);
[err_dte,x] =getErr(p_dte,dte,a,dta,resp,funName);
[err_a_dte,x] =getErr(p_a_dte,dte,a,dta,resp,funName);
[err_a_dta,x] =getErr(p_a_dta,dte,a,dta,resp,funName);
[err_dte_dta,x] =getErr(p_dte_dta,dte,a,dta,resp,funName);
[err_a_dte_dta,x] =getErr(p_a_dte_dta,dte,a,dta,resp,funName);

%% Chi-squared values for just a and dte

% These values replicate R's glm with Anova
% > Anova(glm.out.amp.dte)
% Analysis of Deviance Table (Type II tests)
%
% Response: prob_2
%     LR Chisq Df Pr(>Chisq)
% dte   50.643  1  1.108e-12 ***
% amp    7.872  1   0.005022 **

% degrees of freedom are always 1 because we're always adding 1 factor
chi_a = -2*(err_a_dte-err_dte);
pval_a = 1-chi2cdf(chi_a,1);

chi_dte= -2*(err_a_dte-err_a);
pval_dte = 1-chi2cdf(chi_dte,1);

disp(sprintf('a: chi-squared(%d) = %5.2f, p = %0.7f',...
    1,chi_a,pval_a));

disp(sprintf('dte: chi-squared(%d) = %5.2f, p = %0.7f',...
    1,chi_dte,pval_dte));

%% chi-squared values for a, dte and dta

% These p-values replicate R's glm and 'Anova'
% > Anova(glm.out.amp.dte.dta)
% Analysis of Deviance Table (Type II tests)
%
% Response: prob_2
%     LR Chisq Df Pr(>Chisq)
% amp   7.6109  1   0.005802 **
% dte  18.3959  1  1.794e-05 ***
% dta   7.2029  1   0.007279 **


chi_a = -2*(err_a_dte_dta-err_dte_dta);
pval_a = 1-chi2cdf(chi_a,1);

chi_dte= -2*(err_a_dte_dta-err_a_dta);
pval_dte = 1-chi2cdf(chi_dte,1);

chi_dta= -2*(err_a_dte_dta-err_a_dte);
pval_dta = 1-chi2cdf(chi_dta,1);

disp(sprintf('Including dta:'))

disp(sprintf('a: chi-squared(%d) = %5.2f, p = %0.7f',...
    1,chi_a,pval_a));

disp(sprintf('dte: chi-squared(%d) = %5.2f, p = %0.7f',...
    1,chi_dte,pval_dte));

disp(sprintf('dta: chi-squared(%d) = %5.2f, p = %0.7f',...
    1,chi_dta,pval_dta));

%% Two surfaces - maximal and minimal values of dta.

% Meshgrid for surfaces
surfRez = 21;
binRez = 7;
[aList,dteList] = meshgrid(linspace(50,max(a),surfRez),linspace(min(dte),max(dte),surfRez));

% For fitting, take the best parameters from a_dte and let dta go free

% dta: (for maximal distance to axon)
p_dta = fit('getErr',p_a_dte,{'kdta'},dte,a,dta,resp,funName);
dtaList = dteList;  % set dta to dte (it's maximum)
[~,~,probList_dta] = getErr(p_dta,dteList,aList,dtaList,[],funName);

% dte-dta: (called 'dtema' for minimal distance to axon)
p_dtema = fit('getErr',p_a_dte,{'kdta'},dte,a,dte-dta,resp,funName);
dtaList = dteList;  % again, maximum is dte (where dta = 0)
[~,~,probList_dtema] = getErr(p_dtema,dteList,aList,dtaList,[],funName);

% plot the surfaces
figure(3)
clf
hold on
surf(dteList,aList,probList_dta,'FaceAlpha',.5,'EdgeColor','none');
surf(dteList,aList,probList_dtema,'FaceAlpha',.5,'EdgeColor','none');
set(gca, 'XLim', [0 6000])
set(gca, 'YLim', [0 650])

ylabel('Amplitude (\muA)', 'fontweight','bold');
xlabel('Physical Distance  (\mum)','fontweight','bold');
zlabel('P("2 percepts")', 'fontweight','bold');
grid off
c=  colorbar;
pos  = c.Position;
set(c,'Position',[0.820357142857141 0.19952380952381 0.031499999999999 0.493333333333336])
view(-20,20)

exportgraphics(gcf,['.' filesep 'figures' filesep 'fig7a.eps'],'ContentType','vector')
exportgraphics(gcf,['.' filesep 'figures' filesep 'fig7a.emf'],'ContentType','vector')
%exportgraphics(gcf,['.' filesep 'figures' filesep 'fig7a.pdf'],'ContentType','vector')


%% Contour plots for maximal and minimal dta

contourList = [.65,.75,.85];
f4 = figure(4);
f4.Position = [488 308 385 453];

clf
[c,h]= contour(dteList,aList,probList_dta,contourList);

count = 1;
for i=1:length(contourList)
    n = c(2,count);
    xcdta{i} = c(1,(count+1):(count+n));
    ycdta{i} = c(2,(count+1):(count+n));
    count = count+n+1;
end


[c,h]= contour(dteList,aList,probList_dtema,contourList);

count = 1;
for i=1:length(contourList)
    n = c(2,count);
    xcdtema{i} = c(1,(count+1):(count+n));
    ycdtema{i} = c(2,(count+1):(count+n));
    count = count+n+1;
end

lineStyles = {'g:','g-','g-.'};

clf
hold on
for i=2 %1:length(contourList) % for 75% threshold
    h(i,1) = plot(xcdta{i},ycdta{i},lineStyles{i},'LineWidth',2);
    p = polyfit(xcdta{i},ycdta{i}, 1)
    str{i} = sprintf('%d%%',100*contourList(i));
    text(xcdta{i}(1)-150,ycdta{i}(1)+30,str{i},'HorizontalAlignment','center','Color','g');
    tmp = interp1(ycdta{i}, xcdta{i}, [210 274 476]);
    disp(['predicted spatial res no axon = ', num2str(tmp)])
end

lineStyles = {'b:','b-','b-.'};

for i=2 %1:length(contourList)
    h(i,2) = plot(xcdtema{i},ycdtema{i},lineStyles{i},'LineWidth',2);
    str{i} = sprintf('%d%%',100*contourList(i));
    text(xcdtema{i}(1)+350,ycdtema{i}(1)+30,str{i},'HorizontalAlignment','center','Color','b');
    
end
ylabel('Amplitude (\muA)', 'fontweight','bold');
xlabel('Physical Distance  (\mum)','fontweight','bold');

text(3000,200,'Two Percepts');
text(1000,550,'One Percept');

%find 75% 2p thresh from median thresholds for max axon model
max_axon_val = interp1(yc{2}, xc{2}, [210 274 476]);
disp(['predicted spatial res = ', num2str(tp_val)])
%rectangle('Position',[0 0  tp_val(1) 210],'LineStyle','--', ...
 %         'EdgeColor', [0.5 0.5 0.5 0.5], 'LineWidth',2)
%rectangle('Position',[0 0  tp_val(2) 274],'LineStyle','--', ...
 %         'EdgeColor', [0.5 0.5 0.5 0.5], 'LineWidth',2)
%rectangle('Position',[0 0  tp_val(3) 476],'LineStyle','--', ...
   %       'EdgeColor', [0.5 0.5 0.5 0.5], 'LineWidth',2)

%legend(str,'Location','NorthWest')
set(gca,'XLim',1000*[.575,4])
set(gca,'YLim',[50,700])


exportgraphics(gcf,['.' filesep 'figures' filesep 'fig7b.eps'],'ContentType','vector')
exportgraphics(gcf,['.' filesep 'figures' filesep 'fig7b.emf'],'ContentType','vector')
exportgraphics(gcf,['.' filesep 'figures' filesep 'fig7b.pdf'],'ContentType','vector')
%% Fits, Chi-squared and p-values


grid off

legend(h(2,:),{'Maximal distance to axon','Minimal distance to axon'},...
    'Location','NorthWest');
