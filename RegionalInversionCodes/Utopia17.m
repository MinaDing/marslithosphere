%% now we have the localized data, run MCMC code
clear;
load Utopia_data.mat

Te_range = [0 300];
rhot_range = [2200 3300];
Te_isd = 10;
rhot_isd = 50;
Te_ini = 50;
rhot_ini = 2900;

N = 20000;
burnin = 2000;
gamma = 5;

tic 
%MCMC_ImpactSpatial(Te_range,rhot_range,Te_isd,rhot_isd,Te_ini,rhot_ini,N,burnin,Dvec,FAAvec,FAAvec_sd,la,LM,topox_lm,Grid_d,TOPO_end)
MCMC_ImpactSpatial_coarsened(gamma,Te_range,rhot_range,Te_isd,rhot_isd,Te_ini,rhot_ini,N,burnin,Dvec,FAAvec,FAAvec_sd,la,LM,topox_lm,Grid_d,TOPO_end)
eval(['!mv tmp.mat Utopia_spatial_' num2str(gamma) '.mat']);
toc

%% grid search code


    
%{
    
    %% solve for best-fits 
[rhotrhot, TeTe]=meshgrid(2200:50:3300,0:10:200);
N = numel(rhotrhot);
NRMSE = TeTe*NaN;
for i=1:N
    rhot = rhotrhot(i);
    Te = TeTe(i)*1e3;
    faa = GravityPredict(rhot,Te,Dvec,la,LM,topox_lm,Grid_d,TOPO_end);
    NRMSE(i)=std((FAAvec'-faa)./FAAvec_sd);
end


figure;
[ccc,fff]=contourf(rhotrhot,TeTe,NRMSE,2:0.1:2.5)
clabel(ccc,fff)

hold on
[~,ix]=min(NRMSE);
for i=1:length(ix)
disp(NRMSE(ix(i),i))
plot(rhotrhot(ix(i),i),TeTe(ix(i),i),'ko','Markerfacecolor','k')
end

colormap spring

%}
    
%%

load Utopia_spatial.mat
disp 'Plotting data comparison...'
NRMSE = std((FAAvec'-faa_best)./FAAvec_sd);

string = {['Te=' num2str(Tebest) 'km, sd=' num2str(Tesd) 'km'];
    ['\rhol=' num2str(rhotbest) 'kg/m^3, sd=' num2str(rhotsd) 'kg/m^3'];
    ['NRMSE=', num2str(NRMSE)];
    };

disp(string);

figure;
subplot1(2,1);
subplot1(1)
%plot(Grid_d,Grid_FAA,'k.');
hold on;
%plot(Grid_d,Grid_faa,'r');
errorbar(Dvec,FAAvec,FAAvec_sd,'k')
errorbar(Dvec,faa_best,faa_sd*2,'r')
xlim([0 theta])
text(0,100,string);

subplot1(2)
xlim([0 theta])
errorbar(Dvec,TOPOvec,TOPOvec_sd,'k')
plot(Dvec,fun([7.4 15.9 4.7],Dvec),'k--');

%
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

%

disp 'Plotting covariance after burnin ...'

figure;
set(gcf,'defaultaxescolororder',[0.2 0.2 0.2])

[S,AX,BigAx,H,HAx]=plotmatrix(X(:,burnin:N)','.');
title(BigAx,'Covariance plot of the data')
set(S(:,:),'MarkerSize',3,'MarkerFaceColor',[0.2 0.2 0.2])

set(HAx(1),'Xlim',Te_range)
set(HAx(2),'Xlim',rhot_range)

set(AX(:,1),'Xlim',Te_range,'Xtick',0:50:200)
set(AX(:,2),'Xlim',rhot_range,'Xtick',2500:500:3300)
set(AX(1,:),'Ylim',Te_range,'Ytick',0:50:200)
set(AX(2,:),'Ylim',rhot_range,'Ytick',2500:500:3300)

ylabel(AX(1,1),'Te (km)'); 
ylabel(AX(2,1),'\rhot (kg/m^3)')

xlabel(AX(2,1),'Te (km)');
xlabel(AX(2,2),'\rhot (kg/m^3)')

% now add Gaussion approximation to histogram
axes(HAx(1));
H(1,1).NumBins = 10;
%pd = fitdist(X(1,burnin:N)','Normal');
Dbin = (H(1,1).BinLimits(2)-H(1,1).BinLimits(1))/H(1,1).NumBins;
line(0:200,pdf('normal',0:200,Tebest,Tesd)*(N-burnin+1)*Dbin,'color','b');

axes(HAx(2));
H(1,2).NumBins = 20;
%pd = fitdist(X(2,burnin:N)','Normal');
Dbin = (H(1,2).BinLimits(2)-H(1,2).BinLimits(1))/H(1,2).NumBins;
%line(1000:10:2000,pdf(pd,1000:10:2000)*(N-burnin+1)*Dbin,'color','b');
line(2000:10:3500,pdf('normal',2000:10:3500,rhotbest,rhotsd)*(N-burnin+1)*Dbin,'color','b');

