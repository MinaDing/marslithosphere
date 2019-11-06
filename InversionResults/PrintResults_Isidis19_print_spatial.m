%%
load Isidis_data.mat
load Isidis_spatial_10.mat
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

disp(string);

figure;
subplot1(2,1);
subplot1(1)
%plot(Grid_d,Grid_FAA,'k.');
hold on;
%plot(Grid_d,Grid_faa,'r');
errorbar(Dvec,FAAvec,FAAvec_sd,'k')
errorbar(Dvec,faa_best,faa_sd,'r')
xlim([0 theta])
text(3,100,string);
ylabel 'Free-air anomaly (mGal)'

subplot1(2)
xlim([0 theta])
errorbar(Dvec,TOPOvec,TOPOvec_sd,'k')
plot(Dvec,fun([7.4 15.9 4.7],Dvec),'k--');
set(gca,'xlim',[0 theta],'xtick',0:5:theta,'xticklabel',[0:5:theta]*60)
xlabel 'Distance to basin center (km)'
ylabel 'Topo (km)'

%%
disp 'Plot Marcov Chain...'
figure;
subplot1(2,1);
subplot1(1)
plot(X(1,:),'k');hold on;
ylim(Te_range)
line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
line([burnin N],[Tebest-Tesd Tebest-Tesd],'color','k');
line([burnin N],[Tebest+Tesd Tebest+Tesd],'color','k');
ylabel 'Te (km)'


subplot1(2)
plot(X(2,:),'k');hold on;
ylim(rhot_range)
line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
line([burnin N],[rhotbest-rhotsd rhotbest-rhotsd],'color','k');
line([burnin N],[rhotbest+rhotsd rhotbest+rhotsd],'color','k');
ylabel '\rhot (kg/m^3)'

%%

disp 'Plotting covariance after burnin ...'

subplot1(2,2,'Gap',[0.01 0.01],'YTickL','All');
subplot1(1);
[f,xi,bw] = ksdensity(X(1,burnin:N));
[~,ind]=max(f);
Tebest = xi(ind);
Temin = quantile(X(1,burnin:N),0.16);
Temax = quantile(X(1,burnin:N),0.84);
[counts,centers]=hist(X(1,burnin:N),0:bw:300);
[f,xi,bw] = ksdensity(X(1,burnin:N),0:bw:300);
bar(centers,counts/sum(counts)*100,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
line(xi,f*bw*100,'color','r','linewidth',1);
xlim([0 300])
set(gca,'YMinorTick','on','XMinorTick','on','Xtick',0:50:300,'Ytick',0:2:8,'Layer','top')
ylabel 'Posterior distribution (%)'

subplot1(2);delete(gca);


subplot1(2)
scatter(X(1,burnin:N),X(2,burnin:N),'k.')
xlim([0 300]);
ylim([2200 3300]);
set(gca,'YMinorTick','on','XMinorTick','on','Xtick',0:50:300,'Ytick',2000:250:4000,'Layer','top')
ylabel '\rho_t (kg/m^3)'
xlabel 'Te (km)'

subplot1(3)
[f,xi,bw] = ksdensity(X(2,burnin:N));
[~,ind]=max(f);
rhotbest = xi(ind);
rhotmin = quantile(X(2,burnin:N),0.16);
rhotmax = quantile(X(2,burnin:N),0.84);
[counts,centers]=hist(X(2,burnin:N),2200:bw:3300);
[f,xi,bw] = ksdensity(X(2,burnin:N),2200:bw:3300);
barh(centers,counts/sum(counts),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
line(f*bw,xi,'color','r','linewidth',1);
ylim([2200 3300])
set(gca,'YMinorTick','on','XMinorTick','on','Ytick',2000:250:4000,'Yticklabel',[])
set(gca, 'Layer', 'top')
xlabel 'Posterior distribution (%)'

ezprint('landscape',[])
print(['Isidis19_spatial_cov.pdf'],'-dpdf')
