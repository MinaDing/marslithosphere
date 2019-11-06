clear;clc
capid = 14;

load Volcano14_Elysium_3para_2240_13.mat
load Cap14_admcor.mat;
X(3,:)=-X(3,:);

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

string = {['Te=' num2str(round(Temin)) '-' num2str(round(Temax)) 'km (best = ' num2str(round(Tebest)) 'km)'];
    ['\rhol=' num2str(round(rhotmin)) '-' num2str(round(rhotmax)) 'kg/m^3 (best=' num2str(round(rhotbest)) 'kg/m^3)'];
    ['f=' num2str(round(fmin,2)) '-' num2str(round(fmax,2)) '(best=' num2str(round(fbest,2)) ')'];
    ['reduced chi-square =' num2str(RK2)];
    ['NRMSE=', num2str(NRMSE)];
    };



figure
errorbar(1./DEGS,FADM,FADM_sd,'color',[0.7 0.7 0.7]);hold on;
s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),FADM(s),FADM_sd(s),'k');
hold on;
ylabel('Free-air adm. (mGal/km)')
xlabel 'Spherical harmonic degree'

s = find(DEGS>=degmin&DEGS<=degmax);
errorbar(1./DEGS(s),fadm_best(s),fadm_sd(s),'r')
plot(1./DEGS(s),fadm_best(s),'color','r','linewidth',1);
text(1./(Lwin+5),200,string);

set(gca,'Ytick',0:50:200,'ylim',[25 200],'xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
set(gca,'box','off')

ax2 = axes('Position',ax1_pos);
plot(1./DEGS,FCOR,'b','Parent',ax2);
set(ax2,'Ycolor','b','YAxisLocation','Right','color','none'...
     ,'Ytick',0.5:0.1:1,'Ylim',[0.5 1],'Xtick',1./[70:-10:20],'xlim',[1/(90-Lwin) 1/(Lwin+5)],'xticklabel',[],...
     'YMinorTick','on','XMinorTick','on','xdir','reverse');
set(gca,'box','off')


%%


disp 'Plot Marcov Chain...'

figure;
subplot1(Npara,1);
subplot1(1)
plot(X(1,:),'k');hold on;
line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
line([burnin N],[Tebest-Tesd Tebest-Tesd],'color','k');
line([burnin N],[Tebest+Tesd Tebest+Tesd],'color','k');
ylabel 'Te (km)'

subplot1(2)
plot(X(2,:),'k');hold on;
line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
line([burnin N],[rhotbest-rhotsd rhotbest-rhotsd],'color','k');
line([burnin N],[rhotbest+rhotsd rhotbest+rhotsd],'color','k');
ylabel '\rhot (kg/m^3)'

if Npara==3;
subplot1(3)
plot(X(3,:),'k');hold on;
line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
line([burnin N],[fbest-fsd fbest-fsd],'color','k');
line([burnin N],[fbest+fsd fbest+fsd],'color','k');
ylabel 'f'
end
%%

disp 'Plotting covariance after burnin ...'
Tetick = 0:100:300;
rhotick = 2900:200:3300;
ftick = -1:0.5:1;


figure;
subplot1(3,3,'Gap',[0.01 0.01],'YTickL','All');
subplot1(1);
[f,xi,bw] = ksdensity(X(1,burnin:N));
[~,ind]=max(f);
Tebest = xi(ind);
Temin = quantile(X(1,burnin:N),0.16);
Temax = quantile(X(1,burnin:N),0.84);
[counts,centers]=hist(X(1,burnin:N),Tetick(1):bw:Tetick(end));
[f,xi,bw] = ksdensity(X(1,burnin:N),Tetick(1):bw:Tetick(end));
bar(centers,counts/sum(counts)*100,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
line(xi,f*bw*100,'color','r','linewidth',1);
xlim([Tetick(1) Tetick(end)])
set(gca,'YMinorTick','on','XMinorTick','on','Xtick',Tetick,'Layer','top')
ylabel 'Posterior distribution (%)'

subplot1(2);delete(gca);
subplot1(2);delete(gca);

subplot1(2)
scatter(X(1,burnin:N),X(2,burnin:N),'k.')
xlim([Tetick(1) Tetick(end)]);
ylim([rhotick(1) rhotick(end)]);
set(gca,'YMinorTick','on','XMinorTick','on','Xtick',Tetick,'Ytick',rhotick,'Layer','top')
ylabel '\rho_t (kg/m^3)'

subplot1(3)
[f,xi,bw] = ksdensity(X(2,burnin:N));
[~,ind]=max(f);
rhotbest = xi(ind);
rhotmin = quantile(X(2,burnin:N),0.16);
rhotmax = quantile(X(2,burnin:N),0.84);
[counts,centers]=hist(X(2,burnin:N),rhotick(1):bw:rhotick(end));
[f,xi,bw] = ksdensity(X(2,burnin:N),rhotick(1):bw:rhotick(end));
bar(centers,counts/sum(counts),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
line(xi,f*bw,'color','r','linewidth',1);
xlim([rhotick(1) rhotick(end)])
set(gca,'YMinorTick','on','XMinorTick','on','Xtick',rhotick,'YAxisLocation','Right')
set(gca, 'Layer', 'top')
ylabel 'Posterior distribution (%)'
ylim([0 0.15])

subplot1(4);delete(gca);

subplot1(4)
scatter(X(1,burnin:N),X(3,burnin:N),'k.');
set(gca,'Xtick',Tetick,'XMinorTick','on','Ytick',ftick,'YMinorTick','on','Layer','top')
xlim([Tetick(1) Tetick(end)]);
ylim([ftick(1) ftick(end)]);
ylabel '\alphaf'
xlabel 'Te (km)'

subplot1(5)
scatter(X(2,burnin:N),X(3,burnin:N),'k.');
set(gca,'Xtick',rhotick,'XMinorTick','on','Ytick',ftick,'YMinorTick','on','Layer','top','Yticklabel',[])
xlabel '\rho_t (kg/m^3)'
xlim([rhotick(1) rhotick(end)])
ylim([ftick(1) ftick(end)])

subplot1(6)
[f,xi,bw] = ksdensity(X(3,burnin:N));
[~,ind]=max(f);
rhotbest = xi(ind);
rhotmin = quantile(X(3,burnin:N),0.16);
rhotmax = quantile(X(3,burnin:N),0.84);
[counts,centers]=hist(X(3,burnin:N),ftick(1):bw:ftick(end));
[f,xi,bw] = ksdensity(X(3,burnin:N),ftick(1):bw:ftick(end));
barh(centers,counts/sum(counts),'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
line(f*bw,xi,'color','r','linewidth',1);
ylim([ftick(1) ftick(end)])
set(gca,'YMinorTick','on','XMinorTick','on','Ytick',ftick,'Yticklabel',[])
set(gca, 'Layer', 'top')
xlabel 'Posterir distribution (%)'

ezprint('landscape',[])
print(['/Users/dingmin/Desktop/Cap' num2str(capid) '_cov.pdf'],'-dpdf')

