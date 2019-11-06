function [fadm,fcor,badm,bcor] = localized_synthetic_admcor2(Te,f,alpha,rhol,rhoc,degs,HH,Mij)
% Updated version of localized_synthetic_admcor use Mij for localization
% this function calculate the expected admittance and correlation for
% localized topography HH

% degs should be 0:lmax
% HH is the power spectra of inverted global topography from HH_w 
% Mij is the transformation matrix (lmax-Lwin+1)-by-(lmax+1)
% Lwin is the maximum degree for the localization window
% degs should be 0:lmax

format long
CurrPath = pwd;
addpath([CurrPath '/Subroutines'],'-end');

% Te = 100e3;
% f = 0.5*2900/600;
% alpha = 0.5;
% rhol = 2900;
% degs = 5:100;
 Hc = 50e3;

%%
G = 6.6732e-11;
rhom = 3500;
rhob = rhom-rhoc;
E = 1e11;
nu = 0.25;
Rplanet = 3389.3e3;
gzero = 3.7279;
Robs = 3396e3;


Gt = 4*pi*G*(Rplanet/Robs).^(degs+2).*(1+degs)./(2*degs+1)*1e5;
Gb = 4*pi*G*((Rplanet-Hc)/Robs).^(degs+2).*(1+degs)./(2*degs+1)*1e5;


% if Te == Inf
%     HiHi = degs*0+1;
%     WiWi = f.^2.*HiHi;
%     HiWi = alpha.*f.*HiHi;
%     HW = HiWi;
%     WW = WiWi;
% else
    alphaTOP = 0*degs;
    
    if size(rhol,1)==1
        alphaTOP = alphaT2(degs,rhol,rhom-rhoc,Rplanet,gzero,Te,E,nu);
    else
        for i=1:length(degs)
            alphaTOP(i) = alphaT2(degs(i),rhol(i),rhom-rhoc,Rplanet,gzero,Te,E,nu);
        end

    end
    
    alphaBOT = 1./alphaT2(degs,rhob,rhom-rhob,Rplanet,gzero,Te,E,nu);
    
%     if f == Inf
%         HiHi = degs*0;
%         WiWi = (alphaBOT-1).^2;
%         HiWi = 0;
%         HW = alphaBOT./(alphaBOT-1).^2.*WiWi;
%         WW = alphaBOT.*HW;
%         HHi =degs*0;
%         HiW = degs*0;
%     else

    HiHi = ((1./(1-alphaTOP).^2+f.^2./(alphaBOT-1).^2+2*f.*alpha./(1-alphaTOP)./(alphaBOT-1)).^(-1)).*HH;
    WiWi = f.^2.*HiHi;
    HiWi = alpha.*f.*HiHi;
    HW = alphaTOP./(1-alphaTOP).^2.*HiHi+(alphaTOP+alphaBOT)./(1-alphaTOP)./(alphaBOT-1).*HiWi+alphaBOT./(alphaBOT-1).^2.*WiWi;
    WW = (alphaTOP./(1-alphaTOP)).^2.*HiHi+2*alphaTOP.*alphaBOT./(1-alphaTOP)./(alphaBOT-1).*HiWi+(alphaBOT./(alphaBOT-1)).^2.*WiWi;
    HHi = 1./(1-alphaTOP).*HiHi+1./(alphaBOT-1).*HiWi;
    HiW = alphaTOP./(1-alphaTOP).*HiHi+alphaBOT./(alphaBOT-1).*HiWi;
    

%end

if rhol==rhoc
    
    FH = Gt*rhoc.*HH+Gb.*rhob.*HW;
    FF = rhoc.^2.*Gt.^2.*HH+rhob^2.*Gb.^2.*WW+2*rhoc*rhob*Gt.*Gb.*HW;
    BH = rhob*Gb.*HW;
    BB = rhob^2.*Gb.^2.*WW;
    
else
    FH = Gt*rhoc.*HH+Gt.*(rhol-rhoc).*HHi+Gb.*rhob.*HW;
    FF = rhoc.^2.*Gt.^2.*HH+(rhol-rhoc).^2.*Gt.^2.*HiHi+rhob^2.*Gb.^2.*WW+2*rhoc*(rhol-rhoc).*Gt.^2.*HHi+2*rhoc*rhob*Gt.*Gb.*HW+2*(rhol-rhoc).*rhob.*Gt.*Gb.*HiW;
    BH = (rhol-rhoc).*Gt.*HHi+rhob*Gb.*HW;
    BB = (rhol-rhoc).^2.*Gt.^2.*HiHi+2*(rhol-rhoc)*rhob.*Gt.*Gb.*HiW+rhob^2.*Gb.^2.*WW;
    
end

% localization 

HH_w = Mij*HH;
% or I can use the original HH_w?
FH_w = Mij*FH;
FF_w = Mij*FF;
BH_w = Mij*BH;
BB_w = Mij*BB;


fadm = FH_w./HH_w;
fcor = FH_w./sqrt(FF_w)./sqrt(HH_w);
badm = BH_w./HH_w;
bcor = BH_w./sqrt(BB_w)./sqrt(HH_w);

end