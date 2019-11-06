clear;
capid  = 5;
load RegionalData/AT_admcor.mat

%DEGSi = DEGS;
DEGSi = degs;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load InversionResults/Cap5_AT_3para_2230_10.mat;
degmin = 22;degmax = 30;
disp 'Plotting data comparison...'
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


string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(round(RK2,2))];
    ['NRMS misfit=', num2str(round(NRMSE,2))];
    };
    
figure;
subplot1(2,5,'Gap',[0.007 0.007])
subplot1(1)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
ylabel('Free-air adm. (mGal/km)')
set(gca,'Ytick',0:20:100,'ylim',[0 130],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/(Lwin+7),100,string,'VerticalAlignment','top');
text(0.02,0.98,'a1','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Arabia Terra','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

subplot1(6);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
ylabel('Free-air cor.')
set(gca,'Ytick',0:0.25:1,'Ylim',[0 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'a2','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel 'Spherical harmonic degree'

%%
clear;
capid  = 6;
load RegionalData/NT_admcor.mat

%DEGSi = DEGS;
DEGSi = degs;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load InversionResults/Cap6_NT_3para_2244_20.mat
degmin = 22;degmax = 44;
disp 'Plotting data comparison...'
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


string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(RK2)];
    ['NRMS misfit=', num2str(NRMSE)];
    };
    
subplot1(2)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
set(gca,'Ytick',0:20:100,'ylim',[0 130],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/(Lwin+7),100,string,'VerticalAlignment','top');
text(0.02,0.98,'b1','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Noachis Terra','Units', 'Normalized', 'VerticalAlignment', 'Bottom')



subplot1(7);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
set(gca,'Ytick',0:0.25:1,'Ylim',[0 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'b2','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel 'Spherical harmonic degree'
%%
clear;
capid  = 7;
load RegionalData/TC_admcor.mat

%DEGSi = DEGS;
DEGSi = degs;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load InversionResults/Cap7_TC_3para_4173_25.mat
degmin = 41;degmax = 73;
disp 'Plotting data comparison...'
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


string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(RK2)];
    ['NRMS misfit=', num2str(NRMSE)];
    };
    
subplot1(3)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
set(gca,'Ytick',0:20:100,'ylim',[0 130],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/degmin,100,string,'VerticalAlignment','top');
text(0.02,0.98,'c1','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Terra Cimmeria','Units', 'Normalized', 'VerticalAlignment', 'Bottom')


subplot1(8);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
set(gca,'Ytick',0:0.25:1,'Ylim',[0 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'c2','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel 'Spherical harmonic degree'

%%
load InversionResults/Cap7_TC_3para_2230_10.mat
degmin = 22;degmax = 30;
disp 'Plotting data comparison...'
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


string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(RK2)];
    ['NRMS misfit=', num2str(NRMSE)];
    };
    
subplot1(3)
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
set(gca,'Ytick',0:20:100,'ylim',[0 130],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/degmin,100,string,'VerticalAlignment','top');


subplot1(8);
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
set(gca,'Ytick',0:0.25:1,'Ylim',[0 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])


%%
clear;
capid  = 8;
load RegionalData/TS_admcor.mat

%DEGSi = DEGS;
DEGSi = degs;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load InversionResults/Cap8_Sirenum_3para_2233_15.mat
degmin = 22;degmax = 33;
disp 'Plotting data comparison...'
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


string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(RK2)];
    ['NRMS misfit=', num2str(NRMSE)];
    };
    
subplot1(4)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
set(gca,'Ytick',0:20:100,'ylim',[0 130],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/degmin,100,string,'VerticalAlignment','top');
text(0.02,0.98,'d1','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Terra Sirenum','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

subplot1(9);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
set(gca,'Ytick',0:0.25:1,'Ylim',[0 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'d2','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel 'Spherical harmonic degree'

%%

clear;
capid  = 15;
load RegionalData/Cap15_data.mat

DEGSi = DEGS;
%DEGSi = degs;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load InversionResults/Cap15_SP_3para_15.mat
degmin = 22;degmax = 37;
disp 'Plotting data comparison...'
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


string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) '/' num2str(round(Tebest)) 'km'];
    ['\rho_t=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) '/' num2str(round(rhotbest)) 'kg/m^3'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '/' num2str(round(fbest,2))];
    ['\chi^2/\mu =' num2str(RK2)];
    ['NRMS misfit=', num2str(NRMSE)];
    };
    
subplot1(5)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
set(gca,'Ytick',0:20:100,'ylim',[0 130],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1./(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/degmin,100,string,'VerticalAlignment','top');
text(0.02,0.98,'d1','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Solis Planum','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

subplot1(10);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGSi>=degmin&DEGSi<=degmax);
errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
set(gca,'Ytick',0:0.25:1,'Ylim',[0 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'d2','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel 'Spherical harmonic degree'
%%
%set(gcf,'PaperOrientation','landscape');
ezprint('landscape',[])
%print(['Fig9.pdf'],'-dpdf')