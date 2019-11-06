if 0

% prepare the enviroment 

CurrPath = pwd;
addpath([CurrPath '/Subroutines'],'-end');
format long
load initiation.mat

% add the parameters
capid = 15; 

% load localizaton window coeff. 
% Gridwindow, degwindow = 0:Lwin, Sww
localization = read_sha(strcat('Data/Cap', num2str(capid),'_coef.out'),1);
index = size(localization,1);
[la,ma,~]=ind2lmi(1:index); la=la';
Lwin = max(la);
L = reshape(la',2,index/2)'; L = L(:,1);
MM = reshape(ma,2,index/2)'; MM = MM(:,1);
LM = [L MM];
clear MM L
[Gridwindow,~,~]=plm2xyz([LM reshape(localization',2,index/2)'],1);
Gridwindow = Gridwindow/max(max(Gridwindow));
degwindow = 0:Lwin;
Sww = power_spectra(la,localization,degwindow);

% localize the topography, free-air gravity and Bouguer gravity
[lmcosi,~]=xyz2plm(Grid_topo*1e-3.*Gridwindow,120,'im');
topox_lm = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
index = length(topox_lm);
[la,ma,ia]=ind2lmi(1:index); la=la';
L = reshape(la',2,index/2)'; L = L(:,1);
MM = reshape(ma,2,index/2)'; MM = MM(:,1);
LM = [L MM];
clear MM L
%s = find(la<lmin|la>lmax);topox_lm(s)=0;
%[Grid_topox,~,~]=plm2xyz([LM reshape(topox_lm',2,index/2)'],1);

[lmcosi,~]=xyz2plm(Grid_faa1.*Gridwindow,120,'im'); % use the original gravity field

faax_lm = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
%faax_lm(s)=0;
%[Grid_faax,~,~]=plm2xyz([LM reshape(faax_lm',2,index/2)'],1);

[lmcosi,dw]=xyz2plm(Grid_ba.*Gridwindow,120,'im');
bax_lm = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
%bax_lm(s)=0;
%[Grid_bax,lon,lat]=plm2xyz([LM reshape(bax_lm',2,index/2)'],1);

% now calculate the observed admittance and correlation 
DEGS = Lwin+5:90-Lwin;
[FADM,FCOR,FADM_sd,FCOR_sd]=admcor(la,topox_lm,faax_lm,DEGS);
[BADM,BCOR,BADM_sd,BCOR_sd]=admcor(la,topox_lm,bax_lm,DEGS);

if 0
figure;
subplot1(2,1,'Gap',[0.007 0.007])
subplot1(1);
errorbar(DEGS,FADM,FADM_sd,'k');hold on
ylabel('Free-air adm. (mGal/km)')
set(gca,'Ytick',-100:50:260,'ylim',[0 300],'Xtick',0:20:120,'xlim',[min(DEGS) max(DEGS)],'YMinorTick','on','XMinorTick','on')

subplot1(2);
errorbar(DEGS,FCOR,FCOR_sd,'k');hold on
ylabel('Free-air cor.')
set(gca,'Ytick',-0.5:0.5:1,'Ylim',[-1 1],'Xtick',0:20:120,'xlim',[min(DEGS) max(DEGS)],'YMinorTick','on','XMinorTick','on')
end

% for theoretical curves, assuming Shh = C*l^(-3.5); 
degs = 0:90;
Shh = degs.^(-3.5);
Shh = Shh';
Shh(1)=Shh(2);

% calculate matrix Mij 
Mij = MatrixM(90,Lwin,Sww);

degs = 0:90;
Mij = Mij(Lwin+6:end,:);

Shh = degs.^(-3.5);
Shh = Shh';
Shh(1)=Shh(2);



save Cap15_data.mat


end
%%
load Cap15_data.mat
%FADM_sd2 = sqrt((FADM_sd).^2+5^2);
%FCOR_sd2 = sqrt((FCOR_sd).^2+0.1^2);
FADM_sd2 = FADM_sd;
FCOR_sd2 = FCOR_sd;

%degmin = 27;%Lwin+5;%40
%degmax = 37;%34;%90-Lwin;
degmin = Lwin+5;degmax = 37;
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
rhot_range = [2900 3300];
f_range = [0 10];
alpha_range = [-1 1];

N = 20000; %50,000
burnin = 2000; %5,000

Te_isd = 10;%diff(Te_range)/20;
rhot_isd = 80;%diff(rhot_range)/20;
f_isd = 0.5;%diff(f_range)/20;
alpha_isd = 0.7;%diff(alpha_range)/20;

Te_ini = 120;
rhot_ini = 2990;
f_ini = 3;
alpha_ini = 0;

tic 
%MCMC_flexure_3para3(Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);

gammavec = [1 5 10 15 20 25 30 40];
parfor GI = 1:length(gammavec)
    gamma = gammavec(GI);
    MCMC_flexure_3para3_coarsened(gamma,Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);
end


timeused = toc; 
disp 'MCMC finished'
disp(['Time used is ' num2str(timeused/3600) ' hours']);
%eval(['!move tmp.mat Cap' num2str(capid) '_3para.mat']);
%eval(['!move tmp.mat Cap15_SP_3para_' num2str(degmin) num2str(degmax) '.mat']);
