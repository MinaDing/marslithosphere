CurrPath = pwd;
addpath([CurrPath '/Subroutines'],'-end');
addpath('/Users/dingmin/Documents/MATLAB/Simons');
addpath('/Users/dingmin/Documents/MATLAB/mcgovern/');

clear;
load Isidis_data.mat

[lmcosi,~]=xyz2plm(Grid_topo*1e-3.*Gridwindow,120,'im');
TOPOx_lm = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
index = length(topox_lm);
[la,ma,ia]=ind2lmi(1:index); la=la';
L = reshape(la',2,index/2)'; L = L(:,1);
MM = reshape(ma,2,index/2)'; MM = MM(:,1);
LM = [L MM];
clear MM L

[lmcosi,~]=xyz2plm(Grid_faa1.*Gridwindow,120,'im');
FAAx_lm = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
DEGS = Lwin+5:90-Lwin;
[FADM,FCOR,FADM_sd,FCOR_sd]=admcor(la,TOPOx_lm,FAAx_lm,DEGS);

%%

Te_range = [0 300];
rhot_range = [2200 3300];

Te_isd = 10;
rhot_isd = 50;

Te_ini = 150;
rhot_ini = 2700;

N = 20000;
burnin = 2000;

degmin = 19;
degmax = 27;
s = find(DEGS>=degmin&DEGS<=degmax);

%MCMC_ImpactSpatial(Te_range,rhot_range,Te_isd,rhot_isd,Te_ini,rhot_ini,N,burnin,Dvec,FAAvec,FAAvec_sd,la,LM,topox_lm,Grid_d,TOPO_end);
gammavec = 20;

for i=1:length(gammavec)
    gamma = gammavec(i);
    MCMC_ImpactSpectral_coarsened(gamma,Te_range,rhot_range,Te_isd,rhot_isd,Te_ini,rhot_ini,N,burnin,DEGS(s),FADM(s),FCOR(s),FADM_sd(s),FCOR_sd(s),la,LM,TOPOx_lm,topox_lm,Grid_d,Gridwindow,TOPO_end,Dvec)
    eval(['!mv tmp.mat Isidis19_spectra_' num2str(gamma) '.mat']);
end