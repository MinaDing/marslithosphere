%% now we have the localized data, run MCMC code
clear;
CurrPath = pwd;
addpath('/Users/dingmin/Documents/MATLAB/mcgovern/');
addpath('/Users/dingmin/Documents/MATLAB/Simons');
addpath([CurrPath '/../Subroutines'],'-end');
load Argyre_data.mat

Te_range = [0 300];
rhot_range = [2200 3300];
Te_isd = 10;
rhot_isd = 50;
Te_ini = 100;
rhot_ini = 2900;

N = 20000;
burnin = 2000;

% gammavec = [5 10 15 20 25 30 40 50];
% 
% tic 
% 
% for i=1:length(gammavec)
%     
%     gamma = gammavec(i);
    MCMC_ImpactSpatial(Te_range,rhot_range,Te_isd,rhot_isd,Te_ini,rhot_ini,N,burnin,Dvec,FAAvec,FAAvec_sd,la,LM,topox_lm,Grid_d,TOPO_end)
    eval('!mv tmp.mat Argyre_spatial.mat');
    
    %     MCMC_ImpactSpatial_coarsened(gamma,Te_range,rhot_range,Te_isd,rhot_isd,Te_ini,rhot_ini,N,burnin,Dvec,FAAvec,FAAvec_sd,la,LM,topox_lm,Grid_d,TOPO_end)
%     eval(['!mv tmp.mat Argyre_spatial_' num2str(gamma) '.mat']);
% end

toc
