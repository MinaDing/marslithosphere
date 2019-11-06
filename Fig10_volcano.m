clear;clc
figure
subplot1(3,2,'Gap',[0.007 0.007])

capid = 10;
load RegionalData/Cap10_admcor.mat
load InversionResults/Volcano10_Olympus_3para_2273_30.mat
subplot1(1)

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
set(gca,'Ytick',0:50:200,'ylim',[25 200],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
errorbar(1./DEGS,FCOR,FCOR_sd,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.5:0.1:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],'yticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
set(gca,'box','off')

text(0.08,0.98,string,'Units','Normalized','VerticalAlignment','top');
text(0.02,0.98,'a','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Olympus Mons','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

%%

capid = 9;
load RegionalData/Cap9_admcor.mat
load InversionResults/Volcano9_Alba_3para_2235_15.mat
subplot1(3)

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
set(gca,'Ytick',0:50:200,'ylim',[25 200],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
errorbar(1./DEGS,FCOR,FCOR_sd,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.5:0.1:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],'yticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
set(gca,'box','off')

text(0.08,0.98,string,'Units','Normalized','VerticalAlignment','top');
text(0.02,0.98,'b','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Alba Mons','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

%%
clear;clc
capid = 14;
load RegionalData/Cap14_admcor.mat
load InversionResults/Volcano14_Elysium_3para_2240_13.mat

subplot1(5)

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
set(gca,'Ytick',0:50:200,'ylim',[25 200],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
errorbar(1./DEGS,FCOR,FCOR_sd,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.5:0.1:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],'yticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
set(gca,'box','off')

text(0.08,0.98,string,'Units','Normalized','VerticalAlignment','top');
text(0.02,0.98,'c','Units','Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Elysium Mons','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

%%
clear;clc
capid = 11;
load RegionalData/Cap11_admcor.mat
load InversionResults/Volcano11_Arsia_3para_3164_15.mat

subplot1(2)

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
s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),fadm_best(s),fadm_sd(s),'r')
plot(1./DEGS(s),fadm_best(s),'color','r','linewidth',1);
set(gca,'Ytick',0:50:200,'ylim',[25 200],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20],'yticklabel',[])

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
errorbar(1./DEGS,FCOR,FCOR_sd,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.5:0.1:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
ylabel 'Free-air corr.'
set(gca,'box','off')

text(0.08,0.98,string,'Units','Normalized','VerticalAlignment','top');
text(0.02,0.98,'d','Units','Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Arsia Mons','Units', 'Normalized', 'VerticalAlignment', 'Bottom')


%%
clear;clc
capid = 12;
load RegionalData/Cap12_admcor.mat

subplot1(4)

load InversionResults/Volcano12_Pavonis_3para_3135_4.mat
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
s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),fadm_best(s),fadm_sd(s),'r')
plot(1./DEGS(s),fadm_best(s),'color','r','linewidth',1);
set(gca,'Ytick',0:50:200,'ylim',[25 200],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20],'yticklabel',[])

text(0.08,0.7,string,'Units','Normalized','VerticalAlignment','top');

load InversionResults/Volcano12_Pavonis_3para_4764_14.mat
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

text(0.5,0.7,string,'Units','Normalized','VerticalAlignment','top');

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
errorbar(1./DEGS,FCOR,FCOR_sd,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.5:0.1:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
ylabel 'Free-air corr.'
set(gca,'box','off')

text(0.02,0.98,'e','Units','Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Pavonis Mons','Units', 'Normalized', 'VerticalAlignment', 'Bottom')
%%
clear;clc
capid = 13;
load RegionalData/Cap13_admcor.mat
load InversionResults/Volcano13_Ascraeus_3para_3142_8.mat

subplot1(6)

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
s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),fadm_best(s),fadm_sd(s),'r')
plot(1./DEGS(s),fadm_best(s),'color','r','linewidth',1);
set(gca,'Ytick',0:50:200,'ylim',[25 200],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20],'yticklabel',[])

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
errorbar(1./DEGS,FCOR,FCOR_sd,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.5:0.1:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
ylabel 'Free-air corr.'
xlabel 'Spherical harmonic degree'
set(gca,'box','off')

text(0.08,0.98,string,'Units','Normalized','VerticalAlignment','top');
text(0.02,0.98,'f','Units','Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Ascraeus Mons','Units', 'Normalized', 'VerticalAlignment', 'Bottom')


ezprint('portrait',[])
%print(['Fig11.pdf'],'-dpdf')


