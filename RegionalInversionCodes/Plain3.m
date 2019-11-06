%clear;
% prepare the enviroment 
CurrPath = pwd;
addpath([CurrPath '/Subroutines'],'-end');
format long
load initiation.mat
load VB_admcor.mat
% Shh is the power spectra of the localization window
Sww = Shh;
degwin = 0:Lwin;
DEGS = degs;

% calculate matrix Mij 
Mij = MatrixM(90,Lwin,Sww);

degs = 0:90;
Mij = Mij(Lwin+6:end,:);

Shh = degs.^(-3.5);
Shh = Shh';
Shh(1)=Shh(2);

if 0
figure;
subplot(2,2,1)
errorbar(1./DEGS,FADM,FADM_sd,'k');
ylabel('Free-air adm. (mGal/km)')
%set(gca,'Ytick',-100:50:260,'ylim',[0 100],'Xtick',0:20:120,'xlim',[Lwin+5 max(DEGS)],'YMinorTick','on','XMinorTick','on')
set(gca,'Ytick',-200:100:600,'ylim',[-250 600],'xtick',1./[70:-10:20],'xlim',[min(1./DEGS) 1/22],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
hold on
xlabel 'Spherical harmonic degree'

subplot(2,2,3)
errorbar(1./DEGS,FCOR,FCOR_sd,'k')
ylabel('Free-air cor.')
%set(gca,'Ytick',-0.5:0.5:1,'Ylim',[0.5 1],'Xtick',0:20:120,'xlim',[Lwin+5 max(DEGS)],'YMinorTick','on','XMinorTick','on')
set(gca,'Ytick',0:0.5:1,'Ylim',[-0.5 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGS) 1/22],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
hold on
xlabel 'Spherical harmonic degree'

subplot(2,2,2)
xlabel 'Markov chain iteration'
ylabel 'Te (km)'
ylim([0 200])
hold on;
xlim([0 2000])

subplot(2,2,4)
xlabel 'Markov chain iteraction'
ylabel 'Normalized RMS misfit'
hold on;
xlim([0 2000])
ylim([0 0.7])
end

%%
%FADM_sd2 = sqrt((FADM_sd).^2);
%FCOR_sd2 = sqrt((FCOR_sd).^2);
FADM_sd2 = FADM_sd;
FCOR_sd2 = FCOR_sd;

degmin = Lwin+5;
degmax = 45;
% truncate the observation
s = find(DEGS<degmin|DEGS>degmax);
DEGS(s)=[];
FADM(s)=[];
FCOR(s)=[];
FADM_sd2(s)=[];
FCOR_sd2(s)=[];

[a,b]=size(Mij);
deglocal = 0:a-1;
degs = 0:b-1;


Te_range = [0 300];
rhot_range = [2200 3300];
f_range = [0 20];
alpha_range = [-1 1];

N = 20000; %50,000
burnin = 1000;%500; %5,000


Te_isd = 20;%diff(Te_range)/20;
rhot_isd = 50;%diff(rhot_range)/20;
f_isd = 1;%diff(f_range)/20;
alpha_isd = 0.2;%diff(alpha_range)/20;

Te_ini = 120;
rhot_ini = 3000;
f_ini = 5;
alpha_ini = 0.5;

tic 
%MCMC_flexure_3para3(Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);
%MCMC_flexure_4para(Te_range,rhot_range,f_range,alpha_range,Te_isd,rhot_isd,f_isd,alpha_isd,Te_ini,rhot_ini,f_ini,alpha_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);

% gammavec = [30 40 50 60];
% parfor GI = 1:length(gammavec)
%     disp(GI)
%     gamma = gammavec(GI);
%     MCMC_flexure_4para_coarsened(gamma,Te_range,rhot_range,f_range,alpha_range,Te_isd,rhot_isd,f_isd,alpha_isd,Te_ini,rhot_ini,f_ini,alpha_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);
% end

gammavec = 30;
for GI = 1:length(gammavec)
    gamma = gammavec(GI);
    MCMC_flexure_3para3_coarsened(gamma,Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);
end

timeused = toc; 
disp 'MCMC finished'
disp(['Time used is ' num2str(timeused/3600) ' hours']);
%eval(['!move tmp.mat Cap3_3para.mat']);
%eval(['!move tmp.mat Cap3_VB_4para_' num2str(degmin) num2str(degmax) '.mat']);
