clear;close
load Cap1_data.mat

DEGSi = DEGS;
FADMi = FADM;
FADM_sdi=FADM_sd;
FCORi = FCOR;
FCOR_sdi = FCOR_sd;

load Cap1_NP_3para_FixRhot_15.mat
degmin = 31;
degmax = 48;

% disp 'Plotting data comparison...'
% NRMSE_fadm = std((FADM-fadm_best)./FADM_sd);
% NRMSE_fcor = std((FCOR-fcor_best)./FCOR_sd);

NRMSE = rms([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]);
p = 3;
nu = 2*(degmax-degmin+1)-p;
RK2 = sum(([(FADM-fadm_best)./FADM_sd (FCOR-fcor_best)./FCOR_sd]).^2)/nu;

[f,xi,~] = ksdensity(X(1,burnin:N));
[~,ind]=max(f);
Tebest = xi(ind);
Temin = quantile(X(1,burnin:N),0.16);
Temax = quantile(X(1,burnin:N),0.84);

[f,xi,~] = ksdensity(X(2,burnin:N));
[~,ind]=max(f);
fbest = xi(ind);
fmin = quantile(X(2,burnin:N),0.16);
fmax = quantile(X(2,burnin:N),0.84);

[f,xi,~] = ksdensity(X(3,burnin:N));
[~,ind]=max(f);
alphabest = xi(ind);
alphamin = quantile(X(3,burnin:N),0.16);
alphamax = quantile(X(3,burnin:N),0.84);


string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) 'km (best = ' num2str(round(Tebest)) 'km)'];
    ['f=' num2str(round(fmin)) '-' num2str(round(fmax)) ' (best=' num2str(round(fbest))];
    ['\alpha=' num2str(round(alphamin,2)) '-' num2str(round(alphamax,2)) '(best=' num2str(round(alphabest,2)) ')'];
    ['reduced chi-square =' num2str(RK2)];
    ['NRMSE=', num2str(NRMSE)];
    };


disp(string);
    
figure;
subplot1(2,1)
subplot1(1)
%subplot1(column);
errorbar(1./DEGSi,FADMi,FADM_sdi,'k');hold on;
errorbar(1./DEGS,fadm_best,fadm_sd,'r')
ylabel('Free-air adm. (mGal/km)')
%set(gca,'Ytick',-200:100:600,'ylim',[0 200],'Xtick',0:20:120,'xlim',[Lwin+5 max(DEGSi)],'YMinorTick','on','XMinorTick','on')
set(gca,'Ytick',0:50:100,'ylim',[0 100],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/31],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[])

text(1./(Lwin+5),100,string);

subplot1(2);
errorbar(1./DEGSi,FCORi,FCOR_sdi,'k');hold on;
errorbar(1./DEGS,real(fcor_best),fcor_sd,'r')
ylabel('Free-air cor.')
%set(gca,'Ytick',-0.5:0.5:1,'Ylim',[-0.5 1],'Xtick',0:20:120,'xlim',[Lwin+5 max(DEGSi)],'YMinorTick','on','XMinorTick','on')
set(gca,'Ytick',0:0.5:1,'Ylim',[-0.5 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGSi) 1/31],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])

ezprint('portrait',[])
%print(['Cap' num2str(capid) '.pdf'],'-dpdf')
%%
disp 'Plot Marcov Chain...'


figure;
subplot1(dim,1);
subplot1(1)
plot(X(1,:),'k');hold on;
ylim(Te_range)
line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
line([burnin N],[Tebest-Tesd Tebest-Tesd],'color','k');
line([burnin N],[Tebest+Tesd Tebest+Tesd],'color','k');
ylabel 'Te (km)'


subplot1(2)
plot(X(2,:),'k');hold on;
ylim(f_range)
line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
line([burnin N],[fbest-fsd fbest-fsd],'color','k');
line([burnin N],[fbest+fsd fbest+fsd],'color','k');
ylabel 'f'

subplot1(3)
plot(X(3,:),'k');hold on;
ylim(alpha_range)
line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
line([burnin N],[alphabest-alphasd alphabest-alphasd],'color','k');
line([burnin N],[alphabest+alphasd fbest+alphasd],'color','k');
ylabel 'f'


%%

disp 'Plotting covariance after burnin ...'
figure;
subplot1(3,3,'Gap',[0.01 0.01],'YTickL','All');
subplot1(1);
[f,xi,bw] = ksdensity(X(1,burnin:N));
[~,ind]=max(f);
Tebest = xi(ind);
Temin = quantile(X(1,burnin:N),0.16);
Temax = quantile(X(1,burnin:N),0.84);
[counts,centers]=hist(X(1,burnin:N),0:bw:300);
[f,xi,bw] = ksdensity(X(1,burnin:N),0:bw:300);
bar(centers,counts/sum(counts),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
line(xi,f*bw,'color','r','linewidth',1);
xlim([0 300])
set(gca,'YMinorTick','on','XMinorTick','on','Xtick',0:50:300,'Layer','top')
ylabel 'Posterior distribution (%)'

subplot1(2);delete(gca);
subplot1(2);delete(gca);

subplot1(2)
scatter(X(1,burnin:N),X(2,burnin:N),'k.')
xlim([0 300]);
ylim([0 10]);
set(gca,'YMinorTick','on','XMinorTick','on','Xtick',0:50:300,'Ytick',0:5:10,'Layer','top')
ylabel 'f'

subplot1(3)
[f,xi,bw] = ksdensity(X(2,burnin:N));
[~,ind]=max(f);
fbest = xi(ind);
fmin = quantile(X(2,burnin:N),0.16);
fmax = quantile(X(2,burnin:N),0.84);
[counts,centers]=hist(X(2,burnin:N),0:bw:10);
[f,xi,bw] = ksdensity(X(2,burnin:N),0:bw:10);
bar(centers,counts/sum(counts),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
line(xi,f*bw,'color','r','linewidth',1);
xlim([0 10])
set(gca,'YMinorTick','on','XMinorTick','on','Xtick',0:5:10,'YAxisLocation','Right')
set(gca, 'Layer', 'top')
ylabel 'Posterior distribution (%)'

subplot1(4);delete(gca);

subplot1(4)
scatter(X(1,burnin:N),X(3,burnin:N),'k.');
set(gca,'Xtick',0:100:300,'XMinorTick','on','Ytick',0:0.5:1,'YMinorTick','on','Layer','top')
ylabel '\alpha'
xlabel 'Te (km)'
xlim([0 300])
ylim([0 1])

subplot1(5)
scatter(X(2,burnin:N),X(3,burnin:N),'k.');
set(gca,'Xtick',0:5:10,'XMinorTick','on','Ytick',0:0.5:1,'YMinorTick','on','Layer','top','YTicklabel',[])
xlabel 'f'
xlim([0 10])
ylim([0 1])

subplot1(6)
[f,xi,bw] = ksdensity(X(3,burnin:N));
[~,ind]=max(f);
rhotbest = xi(ind);
rhotmin = quantile(X(3,burnin:N),0.16);
rhotmax = quantile(X(3,burnin:N),0.84);
[counts,centers]=hist(X(3,burnin:N),0:bw:1);
[f,xi,bw] = ksdensity(X(3,burnin:N),0:bw:1);
barh(centers,counts/sum(counts),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
line(f*bw,xi,'color','r','linewidth',1);
ylim([0 1])
set(gca,'YMinorTick','on','XMinorTick','on','Ytick',0:0.5:1,'YTicklabel',[])
set(gca, 'Layer', 'top')
xlabel 'Posterior distribution (%)'

ezprint('landscape',[])
print(['Cap' num2str(capid) '_cov.pdf'],'-dpdf')