load RegionalData/Cap1_data.mat
DEGSi = DEGS;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

%load(['Cap1_NP_3para_3148_15.mat']);
load InversionResults/Cap1_NP_3para_FixRhot_15.mat
degmin = 31;degmax = 48;

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



figure;
subplot1(2,2,'Gap',[0.007 0.007])
subplot1(1)
errorbar(1./DEGSi,FADMi,FADM_sdi,'color',[0.7 0.7 0.7]);hold on;
% s = find(DEGSi>=degmin&DEGSi<=degmax);
% errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
ylabel('Free-air adm. (mGal/km)')
set(gca,'Ytick',0:20:100,'ylim',[0 100],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) max(1./DEGSi)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
%set(gca,'Ytick',-100:50:260,'ylim',[0 100],'Xtick',0:10:120,'xlim',[min(DEGSi) max(DEGSi)],'YMinorTick','on','XMinorTick','on')

text(1/(Lwin+7),100,string,'VerticalAlignment','top');
text(0.02,0.98,'a1','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Northern pole','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

subplot1(3)
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
%s = find(DEGSi>=degmin&DEGSi<=degmax);
%errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
ylabel('Free-air cor.')
%set(gca,'Ytick',-0.5:0.5:1,'Ylim',[0 1],'Xtick',0:10:120,'xlim',[min(DEGSi) max(DEGSi)],'YMinorTick','on','XMinorTick','on')
set(gca,'Ytick',0:0.25:1,'Ylim',[0 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/31],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
xlabel 'Spherical harmonic degree'

text(0.02,0.98,'a2','Units', 'Normalized', 'VerticalAlignment', 'Top')


%%

load('RegionalData/Cap2_data.mat');
DEGSi = DEGS;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

%load(['Cap2_SP_3para_3545_15.mat']);
load InversionResults/Cap2_SP_3para_FixRhot_15.mat
degmin = 35;degmax = 45;

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
% s = find(DEGSi>=degmin&DEGSi<=degmax);
% errorbar(1./DEGSi(s),FADMi(s),FADM_sdi(s),'k');
%errorbar(DEGS,fadm_best,fadm_best-fadm_min',fadm_max'-fadm_best,'b');
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
%ylabel('Free-air adm. (mGal/km)')
%set(gca,'Ytick',-100:50:260,'ylim',[0 100],'Xtick',0:10:120,'xlim',[min(DEGSi) max(DEGSi)],'YMinorTick','on','XMinorTick','on')
set(gca,'Ytick',0:20:100,'ylim',[0 100],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) max(1./DEGSi)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])
text(1/(Lwin+7),100,string,'VerticalAlignment','top');
text(0.02,0.98,'b1','Units', 'Normalized', 'VerticalAlignment', 'Top')
text(0.02,0.02,'Southern pole','Units', 'Normalized', 'VerticalAlignment', 'Bottom')

subplot1(4)
errorbar(1./DEGSi,FCORi,FCOR_sdi,'color',[0.7 0.7 0.7]);hold on;
% s = find(DEGSi>=degmin&DEGSi<=degmax);
% errorbar(1./DEGSi(s),FCORi(s),FCOR_sdi(s),'k');
%errorbar(DEGS,fcor_best,fcor_best-fcor_min',fcor_max'-fcor_best,'b');
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
%ylabel('Free-air cor.')
%set(gca,'Ytick',-0.5:0.5:1,'Ylim',[0 1],'Xtick',0:10:120,'xlim',[min(DEGSi) max(DEGSi)],'YMinorTick','on','XMinorTick','on')
set(gca,'Ytick',0:0.25:1,'Ylim',[0 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/31],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
text(0.02,0.98,'b2','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel 'Spherical harmonic degree'

%ezprint('landscape',[])
%print(['Fig8_2.pdf'],'-dpdf')



