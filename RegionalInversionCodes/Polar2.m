% prepare the enviroment 

if 1
CurrPath = pwd;
addpath([CurrPath '/../Subroutines'])
format long
load initiation.mat

% add the parameters
capid = 2; 

% load localizaton window coeff. 
% Gridwindow, degwindow = 0:Lwin, Sww
localization = read_sha(strcat(CurrPath,'/../RegionalData/Cap', num2str(capid),'_coef.out'),1);
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

%[lmcosi,dw]=xyz2plm(Grid_ba.*Gridwindow,120,'im');
%bax_lm = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
%bax_lm(s)=0;
%[Grid_bax,lon,lat]=plm2xyz([LM reshape(bax_lm',2,index/2)'],1);

% now calculate the observed admittance and correlation 
DEGS = Lwin+5:90-Lwin;
[FADM,FCOR,FADM_sd,FCOR_sd]=admcor(la,topox_lm,faax_lm,DEGS);
%[BADM,BCOR,BADM_sd,BCOR_sd]=admcor(la,topox_lm,bax_lm,DEGS);

if 1
figure;
subplot1(2,1,'Gap',[0.007 0.007])
subplot1(1);
errorbar(DEGS,FADM,FADM_sd,'k');hold on
ylabel('Free-air adm. (mGal/km)')
set(gca,'Ytick',-100:50:260,'ylim',[0 100],'Xtick',0:20:120,'xlim',[min(DEGS) max(DEGS)],'YMinorTick','on','XMinorTick','on')

subplot1(2);
errorbar(DEGS,FCOR,FCOR_sd,'k');hold on
ylabel('Free-air cor.')
set(gca,'Ytick',-0.5:0.5:1,'Ylim',[0 1],'Xtick',0:20:120,'xlim',[min(DEGS) max(DEGS)],'YMinorTick','on','XMinorTick','on')
end

% for theoretical curves, assuming Shh = C*l^(-3.5); 
degs = 0:90;
Shh = degs.^(-3.5);
Shh = Shh';
Shh(1)=Shh(2);

degs_local = Lwin+5:90-Lwin;
Mij = MatrixM(90,Lwin,Sww);
Mij = Mij(Lwin+6:end,:);

% add finite amptitude correction
% [lmcosi,~]=xyz2plm(Grid_FAC.*Gridwindow,120,'im'); % use the original gravity field
% FACx_lm = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
% FAC=admcor(la,topox_lm,FACx_lm,degs_local); %mGal/km/(kg/m^3)

save Cap2_data.mat
end
%%
%addpath('/Users/dingmin/Documents/MATLAB/mcgovern/');
load Cap2_data.mat
%FADM_sd2 = sqrt((FADM_sd).^2+5^2);
%FCOR_sd2 = sqrt((FCOR_sd).^2+0.1^2);
FADM_sd2 = FADM_sd;
FCOR_sd2 = FCOR_sd;

degmin = 35;
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
rhot_range = [1000 2000];
f_range = [0 10];
alpha_range = [-1 1];


Te_isd = 20;%diff(Te_range)/20;
rhot_isd = 50;%diff(rhot_range)/20;
f_isd = 0.1;%diff(f_range)/20;
alpha_isd = 0.1;%diff(alpha_range)/20;

Te_ini = 200;
rhot_ini = 1200;
f_ini = 0.5;
alpha_ini = 0.5;


N = 20000; %50,000
burnin = 1000; %5,000

tic 

gammavec = 15;
for GI = 1:length(gammavec)
    gamma = gammavec(GI);
   % MCMC_flexure_3para3_coarsened(gamma,Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);
    MCMC_flexure_3para_coarsened_FixRhot(gamma,Te_range,f_range,alpha_range,Te_isd,f_isd,alpha_isd,Te_ini,f_ini,alpha_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin);


end



%MCMC_flexure_3para(Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FCOR,FADM_sd2,FCOR_sd2,degs,Shh,Mij,Lwin,FAC);

timeused = toc; 
disp 'MCMC finished'
disp(['Time used is ' num2str(timeused/3600) ' hours']);
%eval(['!move tmp.mat Cap' num2str(capid) '_3para.mat']);
%eval(['!move tmp.mat Cap2_SP_3para_' num2str(degmin) num2str(degmax) '.mat']);
