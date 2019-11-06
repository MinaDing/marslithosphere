# marslithosphere
This repository includes Matlab codes for a journal artical submitted to JGR-Planets: Ding et al., 2019, Variations in Martian Lithospheric Strength Based on Gravity/Topography Analysis

# README #
* This is a repository of the Matlab codes to analyze gravity and topography data of Mars (admittance and correlation analysis) and invert for lithospheric flexural parameters using MCMC method
* Codes to plot the observed and modeled spectra, as well as the Markov chains are also included. 
* For the paper submitted to JGR-Planets: Variations in Martian Lithospheric Strength Based on Gravity/topography Analysis, Ding et al. 

* Please contact Min Ding (dingmin@alum.mit.edu.cn) if you have any question. 
* Created on Jun 14, 2019

### Dependencies ###

* Gravity and topography data: The gravity model JGMRO120d and topography model MarsTopo719 from the Geosciences Node of NASA's Planetary Data System (http://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/) and from the SHTOOLS package (http://sourceforge.net/projects/shtools/), respectively. 

* We used the Matlab codes (alphaB2, alphaT2, cilmn_fac, frigid, gravtabshp4, read_sha, tabload, wiect2p3, gravtabshp4) from Patrick McGovern to model the gravity signature for volcanic montes (Ref: McGovern, P. J., Solomon, S. C., Smith, D. E., Zuber, M. T., Simons, M., Wieczorek, M. A., et al., 2002, Localized gravity/topography admittance and correlation spectra on Mars: Implications for regional and global evolution. JGR-Planets, 19–1. https://doi.org/10.1029/2002JE001854) 

* Matlab codes from Frederik Simons (available online http://geoweb.princeton.edu/people/simons/software.html) are also necessary for data input/output and analysis, including addmon, addmup, ClebschGordan, xyz2plm, plm2xyz, write_sha, Wiegner3j, verbose, shsin, shoos, sha2grid, etc. 

* SHTOOLS (https://shtools.oca.eu/shtools/) is used to find the localization taper for a given cap radius and wavelength.

### Content of This Repository ###

* Fig7_polar.m, Fig8_splain.m, ..., Fig12_Impact.m: a series of Matlab codes to generate Figures 7-12 in the submitted paper. Those figures compare localized admittance/correlation spectra with inverted models using MCMC method. 

* /RegionalData: A series of mat files contain the observed admittance and correlation for all the sub-regions. Called by Fig7_polar.m etc. These mat files were created by the first section of the Matlab codes in /RegionalInversionCodes. SHTOOLS (in fortran or python) is needed to compute localized windows , and sha_read.m is needed to read the spherical harmonic coefficients of the localized windows (presumably the file name is Cap#_coef.out). 

* /RegionalInversionCodes: MCMC inversion codes for each sub-region. Data files in /RegionalData are input for the inversions. 

* /Subroutines: Key functions for the MCMC inversion (see Text S4 of the submitted paper). localized_synthetic_admcor2.m computes the synthetic admittance and correlation functions for plain and polar regions (see Text S1 of the submitted paper). mcgovern_synthetic_admcor.m calls the function gravtabshp4 from Patrick McGovern to forward model the gravity signal for volcanic montes and Valles Marineris. 

* /InversionResults: Results of the MCMC inversion using the codes in /RegionalInversionCodes. Called by Fig7_polar.m etc. Also included are codes (PrintResults_*.m) to print the Markov chains and covariance of the target parameters (in the supplement of the submitted paper). 
