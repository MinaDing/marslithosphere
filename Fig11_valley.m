clear;clc

figure
subplot1(3,1,'Gap',[0.007 0.007])
subplot1(1)

load InversionResults/VallesMarineris1_3para_3251_7.mat
load(['RegionalData/Cap16_admcor1.mat']);

Npara = size(X,1);
NRMSE = rms((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix));
p = 3; 
nu = (degmax-degmin+1)-p;
RK2 = sum(((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix)).^2)/nu;
    
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

[f,xi,bw] = ksdensity(X(3,burnin:N));
[~,ind]=max(f);
fbest = xi(ind);
fmin = quantile(X(3,burnin:N),0.16);
fmax = quantile(X(3,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(round(RK2,2))];
    ['NRMS misfit=', num2str(round(NRMSE,2))];
    };

errorbar(1./DEGS,FADM,FADM_sd,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),FADM(s),FADM_sd(s),'k');
hold on;
ylabel('Free-air adm. (mGal/km)')
xlabel 'Spherical harmonic degree'

s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),fadm_best(s),fadm_sd(s),'r')
plot(1./DEGS(s),fadm_best(s),'color','r','linewidth',1);
set(gca,'Ytick',70:10:110,'ylim',[70 110],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
errorbar(1./DEGS,FCOR,FCOR_sd,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.6:0.2:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
ylabel 'Free-air corr.'
set(gca,'box','off')

text(0.08,0.98,string,'Units','Normalized','VerticalAlignment','top');
text(0.02,0.98,'a','Units', 'Normalized', 'VerticalAlignment', 'Top');
text(0.02,0.02,'Lus Chasma (West)','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

%%
subplot1(2)
load(['RegionalData/Cap16_admcor2.mat']);
load InversionResults/VallesMarineris2_3para_3139_6.mat
Npara = size(X,1);
NRMSE = rms((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix));
p = 3; 
nu = (degmax-degmin+1)-p;
RK2 = sum(((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix)).^2)/nu;
    
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

[f,xi,bw] = ksdensity(X(3,burnin:N));
[~,ind]=max(f);
fbest = xi(ind);
fmin = quantile(X(3,burnin:N),0.16);
fmax = quantile(X(3,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(round(RK2,2))];
    ['NRMS misfit=', num2str(round(NRMSE,2))];
    };

errorbar(1./DEGS,FADM,FADM_sd,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),FADM(s),FADM_sd(s),'k');
hold on;
ylabel('Free-air adm. (mGal/km)')
xlabel 'Spherical harmonic degree'

s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),fadm_best(s),fadm_sd(s),'r')
plot(1./DEGS(s),fadm_best(s),'color','r','linewidth',1);
set(gca,'Ytick',70:10:110,'ylim',[70 110],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.08,0.98,string,'Units','Normalized','VerticalAlignment','top');

load InversionResults/VallesMarineris2_3para_4064_16.mat
Npara = size(X,1);
NRMSE = rms((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix));
p = 3; 
nu = (degmax-degmin+1)-p;
RK2 = sum(((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix)).^2)/nu;
    
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

[f,xi,bw] = ksdensity(X(3,burnin:N));
[~,ind]=max(f);
fbest = xi(ind);
fmin = quantile(X(3,burnin:N),0.16);
fmax = quantile(X(3,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(round(RK2,2))];
    ['NRMS misfit=', num2str(round(NRMSE,2))];
    };

s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),FADM(s),FADM_sd(s),'k');
errorbar(1./DEGS(s),fadm_best(s),fadm_sd(s),'r')
plot(1./DEGS(s),fadm_best(s),'color','r','linewidth',1);
text(0.5,0.98,string,'Units','Normalized','VerticalAlignment','top');

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
errorbar(1./DEGS,FCOR,FCOR_sd,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.6:0.2:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
ylabel 'Free-air corr.'
set(gca,'box','off')


text(0.02,0.98,'b','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Melas Chasma (Middle)','Units', 'Normalized', 'VerticalAlignment', 'Bottom')


%%

subplot1(3)
load(['RegionalData/Cap16_admcor3.mat']);
load InversionResults/VallesMarineris3_3para_3137_7.mat
Npara = size(X,1);
NRMSE = rms((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix));
p = 3; 
nu = (degmax-degmin+1)-p;
RK2 = sum(((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix)).^2)/nu;
    
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

[f,xi,bw] = ksdensity(X(3,burnin:N));
[~,ind]=max(f);
fbest = xi(ind);
fmin = quantile(X(3,burnin:N),0.16);
fmax = quantile(X(3,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(round(RK2,2))];
    ['NRMS misfit=', num2str(round(NRMSE,2))];
    };

errorbar(1./DEGS,FADM,FADM_sd,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),FADM(s),FADM_sd(s),'k');
hold on;
ylabel('Free-air adm. (mGal/km)')
xlabel 'Spherical harmonic degree'

s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),fadm_best(s),fadm_sd(s),'r')
plot(1./DEGS(s),fadm_best(s),'color','r','linewidth',1);
set(gca,'Ytick',70:10:110,'ylim',[70 110],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.08,0.98,string,'Units','Normalized','VerticalAlignment','top');

load InversionResults/VallesMarineris3_3para_3864_18.mat
Npara = size(X,1);
NRMSE = rms((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix));
p = 3; 
nu = (degmax-degmin+1)-p;
RK2 = sum(((FADM(degs_ix)-fadm_best(degs_ix))./FADM_sd(degs_ix)).^2)/nu;
    
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

[f,xi,bw] = ksdensity(X(3,burnin:N));
[~,ind]=max(f);
fbest = xi(ind);
fmin = quantile(X(3,burnin:N),0.16);
fmax = quantile(X(3,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(round(RK2,2))];
    ['NRMS misfit=', num2str(round(NRMSE,2))];
    };

s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),FADM(s),FADM_sd(s),'k');
errorbar(1./DEGS(s),fadm_best(s),fadm_sd(s),'r')
plot(1./DEGS(s),fadm_best(s),'color','r','linewidth',1);
text(0.5,0.98,string,'Units','Normalized','VerticalAlignment','top');

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
errorbar(1./DEGS,FCOR,FCOR_sd,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.6:0.2:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
ylabel 'Free-air corr.'
set(gca,'box','off')


text(0.02,0.98,'c','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Corporates Chasma (East)','Units', 'Normalized', 'VerticalAlignment', 'Bottom')


ezprint('portrait',[])
%print(['Fig12.pdf'],'-dpdf')


