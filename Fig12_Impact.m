% figure 12
% colorp = [255 0 0
% 255 140 0
% 0 255 0
% 0 0 255
% 255 0 255
% 148 0 211]/255;

figure;
subplot1(2,3,'Gap',[0.007 0.007])

%%

load RegionalData/Utopia_data.mat
load InversionResults/Utopia_spatial_5.mat
disp 'Plotting data comparison...'

NRMSE = std((FAAvec'-faa_best)./FAAvec_sd);
p = 3;
nu =length(Dvec)-p;
RK2 = sum(([(FAAvec'-faa_best)./FAAvec_sd]).^2)/nu;

[f,xi,bw] = ksdensity(X(1,burnin:N));
[~,ind]=max(f);
Tebest = xi(ind);
Temin = quantile(X(1,burnin:N),0.16);
Temax = quantile(X(1,burnin:N),0.84);

[f,xi,bw] = ksdensity(X(2,burnin:N));
[~,ind]=max(f);
rhotbest = xi(ind);
rhotmin = quantile(X(2,burnin:N),0.16);
rhotmax = quantile(X(2,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) ' km'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) ' kg/m^3'];
    ['reduced chi-square =' num2str(round(RK2,2))];
    ['NRMSE=', num2str(round(NRMSE,2))];
    };

subplot1(1)
hold on;
errorbar(Dvec,FAAvec,FAAvec_sd,'k')
errorbar(Dvec,faa_best,faa_sd,'r')
ylabel 'Free-air anomaly (mGal)'
set(gca,'ylim',[-200 450],'ytick',-200:200:400,'yminortick','on','xlim',...
    [0 theta],'xtick',0:5:theta,'xticklabel',[],'xminortick','on');

text(0.08,0.5,string,'Units','Normalized','VerticalAlignment','top','color','red');
text(0.02,0.98,'a1','Units', 'Normalized', 'VerticalAlignment', 'Top');
text(0.02,0.02,'Utopia basin','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

subplot1(4)
xlim([0 theta])
errorbar(Dvec,TOPOvec,TOPOvec_sd,'k')
plot(Dvec,fun([7.4 15.9 4.7],Dvec),'k--');
set(gca,'xlim',[0 theta],'xtick',0:5:theta,'xticklabel',[0:5:theta]*60,'ylim',[-8 1],...
    'ytick',-8:2:0,'yminortick','on','xminortick','on');
text(0.02,0.98,'a2','Units', 'Normalized', 'VerticalAlignment', 'Top');
xlabel 'Distance to basin center (km)'
ylabel 'Topo (km)'

%%

load RegionalData/Argyre_data.mat
load InversionResults/Argyre_spatial.mat
disp 'Plotting data comparison...'

NRMSE = std((FAAvec'-faa_best)./FAAvec_sd);
p = 3;
nu =length(Dvec)-p;
RK2 = sum(([(FAAvec'-faa_best)./FAAvec_sd]).^2)/nu;

[f,xi,bw] = ksdensity(X(1,burnin:N));
[~,ind]=max(f);
Tebest = xi(ind);
Temin = quantile(X(1,burnin:N),0.16);
Temax = quantile(X(1,burnin:N),0.84);

[f,xi,bw] = ksdensity(X(2,burnin:N));
[~,ind]=max(f);
rhotbest = xi(ind);
rhotmin = quantile(X(2,burnin:N),0.16);
rhotmax = quantile(X(2,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) ' km'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) ' kg/m^3'];
    ['reduced chi-square =' num2str(round(RK2,2))];
    ['NRMSE=', num2str(round(NRMSE,2))];
    };

subplot1(2)
hold on;
errorbar(Dvec,FAAvec,FAAvec_sd,'k')
errorbar(Dvec,faa_best,faa_sd,'r')
set(gca,'ylim',[-200 450],'ytick',-200:200:400,'yticklabel',[],'yminortick','on','xlim',...
    [0 theta],'xtick',0:5:theta,'xticklabel',[],'xminortick','on');

text(0.08,0.8,string,'Units','Normalized','VerticalAlignment','top','color','red');
text(0.02,0.98,'b1','Units', 'Normalized', 'VerticalAlignment', 'Top');
text(0.02,0.02,'Argyre basin','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

subplot1(5)
hold on
errorbar(Dvec,TOPOvec,TOPOvec_sd,'k')
plot(Dvec,fun([7 7.5 3.75],Dvec),'k--');
set(gca,'xlim',[0 theta],'xtick',0:5:theta,'xticklabel',[0:5:theta]*60,'ylim',[-8 1],...
    'ytick',-8:2:0,'yticklabel',[],'yminortick','on','xminortick','on');
text(0.02,0.98,'b2','Units', 'Normalized', 'VerticalAlignment', 'Top');
xlabel 'Distance to basin center (km)'

%%

load RegionalData/Isidis_data.mat
load InversionResults/Isidis_spatial_10.mat

NRMSE = std((FAAvec'-faa_best)./FAAvec_sd);
p = 3;
nu =length(Dvec)-p;
RK2 = sum(([(FAAvec'-faa_best)./FAAvec_sd]).^2)/nu;

[f,xi,bw] = ksdensity(X(1,burnin:N));
[~,ind]=max(f);
Tebest = xi(ind);
Temin = quantile(X(1,burnin:N),0.16);
Temax = quantile(X(1,burnin:N),0.84);

[f,xi,bw] = ksdensity(X(2,burnin:N));
[~,ind]=max(f);
rhotbest = xi(ind);
rhotmin = quantile(X(2,burnin:N),0.16);
rhotmax = quantile(X(2,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) ' km'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) ' kg/m^3'];
    ['reduced chi-square =' num2str(round(RK2,2))];
    ['NRMSE=', num2str(round(NRMSE,2))];
    };

subplot1(3)
hold on;
errorbar(Dvec,FAAvec,FAAvec_sd,'k')
errorbar(Dvec,faa_best,faa_sd,'r')
set(gca,'ylim',[-200 450],'ytick',-200:200:400,'yticklabel',[],'yminortick','on','xlim',...
    [0 theta],'xtick',0:5:theta,'xticklabel',[],'xminortick','on');

text(0.08,0.8,string,'Units','Normalized','VerticalAlignment','top','color','red');
text(0.02,0.98,'c1','Units', 'Normalized', 'VerticalAlignment', 'Top');
text(0.02,0.02,'Isidis basin','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

subplot1(6)
hold on
errorbar(Dvec,TOPOvec,TOPOvec_sd,'k')
plot(Dvec,fun([7 7.5 3.75],Dvec),'k--');
set(gca,'xlim',[0 theta],'xtick',0:5:theta,'xticklabel',[0:5:theta]*60,'ylim',[-8 1],...
    'ytick',-8:2:0,'yticklabel',[],'yminortick','on','xminortick','on');
text(0.02,0.98,'c2','Units', 'Normalized', 'VerticalAlignment', 'Top');
xlabel 'Distance to basin center (km)'

ezprint('landscape',[]);
%print('Fig12_impact.pdf','-dpdf')

