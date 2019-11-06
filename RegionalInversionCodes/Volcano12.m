% to use the new PC to run this model, I verified this to solve for a
% specific region
if 0
format long
capid = 12; 
rhobc = 2900;
Rplanet = 3389.5e3; 	%[m]
robs = 3396e3;
lmax = 120;
dlonlat = 1;
G   = 6.6732e-11;		%N m^2 / Kg^2, universal gravity constant
GM = 0.4282837452600000e14;
M = GM/G;
gzero = GM/(Rplanet.^2); 	%[m/s^2] used for loading calculations
Tc = 50e3;
zb = Tc;
rhoc = 2900;
rhom = 3500;
rhob = rhom-rhoc; 
E = 1.0e11;         % [Pa]
nu = 0.25;          % [ ]

% load observed spectra
% load spherical harmonic shape and geoid expansion in Oded's format
% Silm: shape file with planetary radius.
% Hilm: 'topo' file with geoid-referenced topography.
% watch out for difficulties with differing reference radii for shape and 
% gravity expansions

Silm = read_sha(['Data/MarsTopo120_noheader.shape'],1);
index=size(Silm,1);
[la,ma,~]=ind2lmi(1:index); la=la';
L = reshape(la',2,index/2)'; L = L(:,1);
MM = reshape(ma,2,index/2)'; MM = MM(:,1);
LM = [L MM];
clear MM L ma 
%Silm(7:8) = Silm(7:8)*0.05;
Gpilm = read_sha(['Data/jgmro_120d_sha_noheader.tab'],1);
%Gpilm(7:8) = 0.05*Gpilm(7:8);
Geoilm = Gpilm*3396e3; 
Geoilm(1) = Silm(1);
Hilm = Silm-Geoilm;
% supply expansion with only 95 percent of J2 (hydrostatic contribution)
% this contribution will be added back to the deflected surfaces so that 
%we don't model the hydrostatic flattening as a load.
hydrofac = 0.95;
Geohydroilm = zeros(size(Geoilm));
Geohydroilm(7) = hydrofac.*Geoilm(7); % entry 7 is the J2 term
% set nmax = number of terms in Wieczorek and Phillips expansion
% By the way, nmax = 1 is the mass sheet approximation.
%nmax = 10;
nmax = 10; % this is fine for degrees up to 90 or so.
% initialize the cofln matrix for the chosen size of the harmonic expansion 
% by calling "cilmn_fac.m"
% this matrix consists of n column vectors
for ic = 1:nmax
    cofln(:,ic) = cilmn_fac(Hilm,ic);
end
[~,cilmpa] = wiect2p3(Silm,Rplanet,lmax,dlonlat,G,M,rhobc,nmax,cofln);
upconsurf = (Rplanet/robs).^la;
BCilm=cilmpa.*upconsurf.*(la+1).*1e5*G*M/robs^2;
FAAilm=Gpilm.*upconsurf.*(la+1).*1e5*G*M/robs^2;
BAilm = FAAilm-BCilm;

% localization
% load window coefficients
localization = read_sha(strcat('Data/Cap', num2str(capid),'_coef2.out'),1);
index2 = size(localization,1);
[la2,ma2,~]=ind2lmi(1:index2); la2=la2';
L2 = reshape(la2',2,index2/2)'; L2 = L2(:,1);
MM2 = reshape(ma2,2,index2/2)'; MM2 = MM2(:,1);
LM2 = [L2 MM2];
Lwin = max(la2);
clear MM2 L2 la2 ma2
[Gridwindow,~,~]=plm2xyz([LM2 reshape(localization',2,index2/2)'],1);
Gridwindow = Gridwindow/max(max(Gridwindow));
% degwindow = 0:Lwin;
% Shh = power_spectra(la,localization,degwindow);

% localize the topography
[Grid_topo,~,~]=plm2xyz([LM reshape(Hilm',2,index/2)'],1);
[lmcosi,~]=xyz2plm(Grid_topo.*Gridwindow,120,'im');
Hilm_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';

% localize the observed gravity 
[Grid_FAA,~,~]=plm2xyz([LM reshape(FAAilm',2,index/2)'],1);
[lmcosi,~]=xyz2plm(Grid_FAA.*Gridwindow,120,'im');
FAAilm_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';

[Grid_BA,~,~]=plm2xyz([LM reshape(BAilm',2,index/2)'],1);
[lmcosi,~]=xyz2plm(Grid_BA.*Gridwindow,120,'im');
BAilm_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';

degs = Lwin+1:90-Lwin;

[FADM,FCOR,FADM_sd,FCOR_sd]=admcor(la,Hilm_w,FAAilm_w,degs);
FADM = FADM*1e3;FADM_sd = FADM_sd*1e3;
[~,BCOR,~,BCOR_sd]=admcor(la,Hilm_w,BAilm_w,degs);

save(['Cap' num2str(capid) '_admcor.mat'],'degs','FADM','FADM_sd','FCOR','FCOR_sd','BCOR','BCOR_sd','Lwin','la','LM','BCilm','Grid_topo','Gridwindow','Hilm','Silm','Geohydroilm','Rplanet','robs','lmax','dlonlat','M','G','gzero','Tc','rhoc','rhom','rhob','zb','E','nu','cofln','nmax');


figure;
subplot1(3,1,'Gap',[0.007 0.007])
subplot1(1);
errorbar(degs,FADM,FADM_sd,'k');hold on
ylabel('Free-air adm. (mGal/km)')
set(gca,'Ytick',-100:50:260,'ylim',[100 200],'Xtick',0:20:120,'xlim',[Lwin+5 max(degs)],'YMinorTick','on','XMinorTick','on')

subplot1(2);
errorbar(degs,FCOR,FCOR_sd,'k');hold on
ylabel('Free-air cor.')
set(gca,'Ytick',-0.5:0.5:1,'Ylim',[0.8 1],'Xtick',0:20:120,'xlim',[Lwin+5 max(degs)],'YMinorTick','on','XMinorTick','on')

end

% %%
% 
% if 0
% clear;clc
% disp 'now test the situation for surface load only'
% maxNumCompThreads(1);
% capid = 9;
% load(['Cap' num2str(capid) '_admcor.mat']);
% a = 2;
% b = 5; % 10 mGal/km
% FADM_sd2 = sqrt((FADM_sd*a).^2+b^2);
% 
% degmin = 20;
% degmax = 35;
% 
% tic                                                    
% MCMC_flexure_VolcanicMontes_tponly([0 200],[2900 3200],10,25,2000,200,degs,FADM,FADM_sd2,capid,degmin,degmax);
% timeused = toc; 
% disp 'MCMC finished'
% disp(['Time used is ' num2str(timeused/3600) ' hours']);
% !move tmp.mat Cap9_tponly.mat
% 
% end

%% now we have the localized data, run MCMC code for positive f situation

%clear;

capid = 12;
load(['Cap' num2str(capid) '_admcor.mat']);
FADM_sd2 = FADM_sd;

Te_range = [0 300];
rhot_range = [2900 3300];
f_range = [-1 1];

Te_isd = 10;
rhot_isd = 25;
f_isd = 0.08;

Te_ini = 20;
rhot_ini = 3200;
f_ini = 0.1;

N = 10000;
burnin = 1000;

tic


% degmin = 31; degmax = 35;
% gammavec = [4];
% for GI = 1:length(gammavec)
%     gamma = gammavec(GI);
%     MCMC_flexure_VolcanicMontes_coarsened(gamma,Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,degs,FADM,FADM_sd2,capid,degmin,degmax);
% end

degmin = 47;degmax = 64;
gammavec = [14];
for GI = 1:length(gammavec)
    gamma = gammavec(GI);
    MCMC_flexure_VolcanicMontes_coarsened(gamma,Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,degs,FADM,FADM_sd2,capid,degmin,degmax);
end


timeused = toc; 
disp 'MCMC finished'
disp(['Time used is ' num2str(timeused/60) ' minutes']);


%MCMC_flexure_VolcanicMontes(Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,degs,FADM,FADM_sd2,capid,degmin,degmax);

%eval(['!move tmp.mat Cap' num2str(capid) '.mat']);
%eval(['!move tmp.mat Volcano12_Pavonis_3para_' num2str(degmin) num2str(degmax) '.mat']);

% %% negative f situation
% clear;
% capid = 9;
% degmin = 20;
% degmax = 35;
% a = 2;
% b = 10; % mGal/km
% 
% tic
% load(['Cap' num2str(capid) '_admcor.mat']);
% FADM_sd2 = sqrt((FADM_sd*a).^2+b^2); 
% MCMC_flexure_VolcanicMontes([0 100],[2900 3200],[-0.3 0],10,10,0.03,3000,400,degs,FADM,FADM_sd2,capid,degmin,degmax);
% timeused = toc; 
% disp 'MCMC finished'
% disp(['Time used is ' num2str(timeused/3600) ' hours']);
% eval(['!move tmp.mat Cap' num2str(capid) '_negative.mat']);


