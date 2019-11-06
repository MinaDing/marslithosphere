clear;

figure;
subplot1(4,2,'Gap',[0.007 0.007])

%%
capid  = 3;
load RegionalData/VB_admcor.mat

DEGSi = degs;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load InversionResults/Cap3_VB_4para_2245_40.mat
degmin = 22;degmax = 45;
NRMSE = rms([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]);
p = 4;
nu = 2*(degmax-degmin+1)-p;
RK2 = sum(([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]).^2)/nu;

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

[f,xi,bw] = ksdensity(X(4,burnin:N));
[~,ind]=max(f);
alphabest = xi(ind);
alphamin = quantile(X(4,burnin:N),0.16);
alphamax = quantile(X(4,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) ' km'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) ' kg/m^3'];
    ['f=' num2str(round(fmin,1)) '-' num2str(round(fmax,1)) '/' num2str(round(fbest,1))];
    ['\alpha=' num2str(round(alphamin,1)) '-' num2str(round(alphamax,1)) '/' num2str(round(alphabest,1))];    
    ['reduced chi-square =' num2str(round(RK2,2))];
    ['NRMSE=', num2str(round(NRMSE,2))];
    };
    

subplot1(1+4)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
ylabel('Free-air adm. (mGal/km)')
set(gca,'Ytick',-200:200:600,'ylim',[-250 600],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/(Lwin+7),400,string,'VerticalAlignment','top');
text(0.02,0.98,'a3','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot1(3+4);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
ylabel('Free-air cor.')
set(gca,'Ytick',-0.5:0.5:1,'Ylim',[-0.5 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'a4','Units', 'Normalized', 'VerticalAlignment', 'Top')

%%
clear;
capid  = 4;
load RegionalData/AP_admcor.mat

%DEGSi = DEGS;
DEGSi = degs;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load(['InversionResults/Cap4_AP_4para_2228_15.mat']);
degmin = 22; degmax = 28;
NRMSE = rms([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]);
p = 4;
nu = 2*(degmax-degmin+1)-p;
RK2 = sum(([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]).^2)/nu;

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

[f,xi,bw] = ksdensity(X(4,burnin:N));
[~,ind]=max(f);
alphabest = xi(ind);
alphamin = quantile(X(4,burnin:N),0.16);
alphamax = quantile(X(4,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) ' km'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) ' kg/m^3'];
    ['f=' num2str(round(fmin,1)) '-' num2str(round(fmax,1)) '/' num2str(round(fbest,1))];
    ['\alpha=' num2str(round(alphamin,1)) '-' num2str(round(alphamax,1)) '/' num2str(round(alphabest,1))];    
    ['reduced chi-square =' num2str(round(RK2,2))];
    ['NRMSE=', num2str(round(NRMSE,2))];
    };

subplot1(2+4)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
set(gca,'Ytick',-200:200:600,'ylim',[-250 600],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/degmin,400,string,'VerticalAlignment','top');
text(0.02,0.98,'b3','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot1(4+4);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
set(gca,'Ytick',-0.5:0.5:1,'Ylim',[-0.5 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'b4','Units', 'Normalized', 'VerticalAlignment', 'Top')

%%
load(['InversionResults/Cap4_AP_4para_3751_15.mat']);
degmin = 37; degmax = 51;
NRMSE = rms([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]);
p = 4;
nu = 2*(degmax-degmin+1)-p;
RK2 = sum(([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]).^2)/nu;

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

[f,xi,bw] = ksdensity(X(4,burnin:N));
[~,ind]=max(f);
alphabest = xi(ind);
alphamin = quantile(X(4,burnin:N),0.16);
alphamax = quantile(X(4,burnin:N),0.84);

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) ' km'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) ' kg/m^3'];
    ['f=' num2str(round(fmin,1)) '-' num2str(round(fmax,1)) '/' num2str(round(fbest,1))];
    ['\alpha=' num2str(round(alphamin,1)) '-' num2str(round(alphamax,1)) '/' num2str(round(alphabest,1))];    
    ['reduced chi-square =' num2str(round(RK2,2))];
    ['NRMSE=', num2str(round(NRMSE,2))];
    };

subplot1(2+4)
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
text(1/degmin,400,string,'VerticalAlignment','top');


subplot1(4+4);
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')

%%
capid  = 3;
load RegionalData/VB_admcor.mat

DEGSi = degs;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load InversionResults/Cap3_VB_3para_2245_30.mat
degmin = 22;degmax = 45;
NRMSE = rms([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]);
p = 3;
nu = 2*(degmax-degmin+1)-p;
RK2 = sum(([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]).^2)/nu;

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

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) ' km'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) ' kg/m^3'];
    ['f=' num2str(round(fmin,1)) '-' num2str(round(fmax,1)) '/' num2str(round(fbest,1))];
    ['reduced chi-square =' num2str(round(RK2,2))];
    ['NRMSE=', num2str(round(NRMSE,2))];
    };
    

subplot1(1)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
ylabel('Free-air adm. (mGal/km)')
set(gca,'Ytick',-200:200:600,'ylim',[-250 600],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/(Lwin+7),400,string,'VerticalAlignment','top');
text(0.02,0.98,'a1','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot1(3);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
ylabel('Free-air cor.')
set(gca,'Ytick',-0.5:0.5:1,'Ylim',[-0.5 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'a2','Units', 'Normalized', 'VerticalAlignment', 'Top')


%%
clear;
capid  = 4;
load RegionalData/AP_admcor.mat

%DEGSi = DEGS;
DEGSi = degs;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load InversionResults/Cap4_AP_3para_2228_15.mat
degmin = 22; degmax = 28;
NRMSE = rms([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]);
p = 3;
nu = 2*(degmax-degmin+1)-p;
RK2 = sum(([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]).^2)/nu;

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


string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) ' km'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) ' kg/m^3'];
    ['f=' num2str(round(fmin,1)) '-' num2str(round(fmax,1)) '/' num2str(round(fbest,1))];
    ['reduced chi-square =' num2str(round(RK2,2))];
    ['NRMSE=', num2str(round(NRMSE,2))];
    };
    
subplot1(2)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
set(gca,'Ytick',-200:200:600,'ylim',[-250 600],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/degmin,400,string,'VerticalAlignment','top');
text(0.02,0.98,'b3','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot1(4);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
set(gca,'Ytick',-0.5:0.5:1,'Ylim',[-0.5 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'b4','Units', 'Normalized', 'VerticalAlignment', 'Top')

%%
load InversionResults/Cap4_AP_3para_3751_25.mat
degmin = 37; degmax = 51;
NRMSE = rms([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]);
p = 3;
nu = 2*(degmax-degmin+1)-p;
RK2 = sum(([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]).^2)/nu;

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

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) ' km'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) ' kg/m^3'];
    ['f=' num2str(round(fmin,1)) '-' num2str(round(fmax,1)) '/' num2str(round(fbest,1))];
    ['reduced chi-square =' num2str(round(RK2,2))];
    ['NRMSE=', num2str(round(NRMSE,2))];
    };

subplot1(2)
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
text(1/degmin,400,string,'VerticalAlignment','top');


subplot1(4);
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')


%%
ezprint('portrait',[])
%set(gcf,'PaperOrientation','landscape');
%print('Fig10.pdf','-dpdf')