% prepare the enviroment 
CurrPath = pwd;
addpath([CurrPath '/Subroutines'],'-end');
format long
load initiation.mat

capid = 7;
load TC_admcor.mat
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
subplot1(2,1)
subplot1(1);
errorbar(DEGS,FADM,FADM_sd,'k');
ylabel('Free-air adm. (mGal/km)')
set(gca,'Ytick',-100:50:260,'ylim',[0 100],'Xtick',0:20:120,'xlim',[Lwin+5 max(DEGS)],'YMinorTick','on','XMinorTick','on')

subplot1(2);
errorbar(DEGS,FCOR,FCOR_sd,'k')
ylabel('Free-air cor.')
set(gca,'Ytick',-0.5:0.5:1,'Ylim',[0.5 1],'Xtick',0:20:120,'xlim',[Lwin+5 max(DEGS)],'YMinorTick','on','XMinorTick','on')
end

%%
%FADM_sd2 = sqrt((FADM_sd).^2+5^2);
%FCOR_sd2 = sqrt((FCOR_sd).^2+0.05^2);
FADM_sd2 = FADM_sd;
FCOR_sd2 = max(FCOR_sd,0.05);

%degmin = Lwin+5;degmax = 32;
degmin = 41; degmax = 90-Lwin;

%degmin = 41;%Lwin+5;%40
%degmax = 70;%34;%90-Lwin;
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
f_range = [0 10];
alpha_range = [-1 1];

N = 20000; %50,000
burnin = 2000; %5,000

Te_isd = 10;%diff(Te_range)/20;
rhot_isd = 100;%diff(rhot_range)/20;
f_isd = 0.5;%diff(f_range)/20;
alpha_isd = 0.1;%diff(alpha_range)/20;

Te_ini = 50;
rhot_ini = 2500;
f_ini = 1;
alpha_ini = 0;

tic 

gammavec = [1 5 10 15 20 25 30 40];
parfor GI = 1:length(gammavec)
    gamma = gammavec(GI);
    MCMC_flexure_3para3_coarsened(gamma,Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);
end

%MCMC_flexure_3para3(Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);

timeused = toc; 
disp 'MCMC finished'
disp(['Time used is ' num2str(timeused/3600) ' hours']);
%eval(['!move tmp.mat Cap' num2str(capid) '_3para.mat']);
%eval(['!move tmp.mat Cap7_TC_3para_' num2str(degmin) num2str(degmax) '.mat']);