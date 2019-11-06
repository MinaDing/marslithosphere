%clear;
% prepare the enviroment 
CurrPath = pwd;
addpath([CurrPath '/Subroutines'],'-end');
format long
load initiation.mat
load AP_admcor.mat
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

if 1
figure;
subplot(2,1,1)
errorbar(DEGS,FADM,FADM_sd,'k');
ylabel('Free-air adm. (mGal/km)')
%set(gca,'Ytick',-100:50:260,'ylim',[0 100],'Xtick',0:20:120,'xlim',[Lwin+5 max(DEGS)],'YMinorTick','on','XMinorTick','on')
%set(gca,'Ytick',-200:100:600,'ylim',[-250 600],'xtick',1./[70:-10:20],'xlim',[min(1./DEGS) 1/22],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
hold on
xlabel 'Spherical harmonic degree'

subplot(2,1,2)
errorbar(DEGS,FCOR,FCOR_sd,'k')
ylabel('Free-air cor.')
%set(gca,'Ytick',-0.5:0.5:1,'Ylim',[0.5 1],'Xtick',0:20:120,'xlim',[Lwin+5 max(DEGS)],'YMinorTick','on','XMinorTick','on')
%set(gca,'Ytick',0:0.5:1,'Ylim',[-0.5 1],'xtick',1./[70:-10:20],'xlim',[min(1./DEGS) 1/22],'YMinorTick','on','XMinorTick','on','xdir','reverse','xticklabel',[70:-10:20])
hold on
xlabel 'Spherical harmonic degree'

end

%%


FADM_sd2 = FADM_sd;
FCOR_sd2 = FCOR_sd;


[a,b]=size(Mij);
deglocal = 0:a-1;
degs = 0:b-1;


Te_range = [0 300];
rhot_range = [2200 3300];
f_range = [0 15];
alpha_range = [-1 1];

N = 20000; %50,000
burnin = 2000;%1000; %5,000

Te_isd = 10;%diff(Te_range)/20;
rhot_isd = 200;%diff(rhot_range)/20;
f_isd = 0.5;%diff(f_range)/20;
alpha_isd = 0.1;%diff(alpha_range)/20;

Te_ini = 100;
rhot_ini = 3000;
f_ini = 4;
alpha_ini = -0.2;

tic 
%MCMC_flexure_3para3(Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);
%MCMC_flexure_4para(Te_range,rhot_range,f_range,alpha_range,Te_isd,rhot_isd,f_isd,alpha_isd,Te_ini,rhot_ini,f_ini,alpha_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);


gammavec = [25];

%degmin = 22;degmax = 28;
degmin = 37;degmax = 51;
s = find(DEGS<degmin|DEGS>degmax);
DEGS(s)=[];
FADM(s)=[];
FCOR(s)=[];
FADM_sd2(s)=[];
FCOR_sd2(s)=[];

for GI = 1:length(gammavec)
    gamma = gammavec(GI);
    %MCMC_flexure_4para_coarsened(gamma,Te_range,rhot_range,f_range,alpha_range,Te_isd,rhot_isd,f_isd,alpha_isd,Te_ini,rhot_ini,f_ini,alpha_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);
    MCMC_flexure_3para3_coarsened(gamma,Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);

end


timeused = toc; 
disp 'MCMC finished'
disp(['Time used is ' num2str(timeused/3600) ' hours']);
%eval(['!move tmp.mat Cap4_3para.mat']);
%eval(['!move tmp.mat Cap4_AP_4para_' num2str(degmin) num2str(degmax) '.mat']);
