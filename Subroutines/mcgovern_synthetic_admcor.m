function [fadm,fcor] = mcgovern_synthetic_admcor(Te,rhol,f,capid,degrees)
% variables
defval('capid',10);
defval('Te',100e3);% meters
defval('rhol',3000); % kg/m^3 
defval('f',-0.2);
defval('degrees',18:72);
load(['Cap' num2str(capid) '_admcor.mat']);

%delrhomoho = rhom - rhoc;  %[kg/m^3] density contrast at moho (for gravity calc)
%delrhosurf = rhoc - rhol;  %[kg/m^3] density contrast between load and crust  
				%(for gravity calc)
rhorestore = rhom - rhol;  %[kg/m^3] sub-lithosphere restoring force 
				% = rhom if no infill assumed
				% = rhom - rhol if load infill assumed 

%% Model the gravity signals
% call script that produces the gravity and potential coefficients 
[~,cilmpall,~,~,Ht,Wt,Hb,Wb] = gravtabshp4(Hilm,Silm,Geohydroilm,Rplanet,robs,lmax,dlonlat,M,G,gzero,Te,Tc,rhol,rhoc,rhom,rhorestore,rhob,zb,f,E,nu,cofln,nmax);
upconsurf = (Rplanet/robs).^la;
faailm=cilmpall.*upconsurf.*(la+1).*1e5*G*M/robs^2;
bailm = faailm-BCilm;
index = length(faailm);
% apply localization 

[Grid_faa,~,~]=plm2xyz([LM reshape(faailm',2,index/2)'],1);
[lmcosi,~]=xyz2plm(Grid_faa.*Gridwindow,120,'im');
faailm_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
% 
% [Grid_ba,~,~]=plm2xyz([LM reshape(bailm',2,index/2)'],1);
% [lmcosi,~]=xyz2plm(Grid_ba.*Gridwindow,120,'im');
% bailm_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';

[lmcosi,~]=xyz2plm(Grid_topo.*Gridwindow,120,'im');
Hilm_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';

[fadm,fcor,~,~]=admcor(la,Hilm_w,faailm_w,degrees);
fadm = fadm*1e3;
% [badm,bcor,~,~]=admcor(la,Hilm_w,bailm_w,degrees);
% badm = badm*1e3;


%%

if 0
    figure
    subplot1(3,1);
    subplot1(1)
    hold on;plot(degrees,fadm);
    subplot1(2)
    hold on;plot(degrees,fcor);
    subplot1(3);
    hold on; plot(degrees,bcor);

%%

if 0
lon = 0:360;
lat = 90:-1:-90;
[lonlon,latlat]=meshgrid(lon,lat);  
[Grid_Ht,~,~]=plm2xyz([LM reshape(Ht',2,index/2)'],1);
[Grid_Wt,~,~]=plm2xyz([LM reshape(Wt',2,index/2)'],1);
[Grid_Hb,~,~]=plm2xyz([LM reshape(Hb',2,index/2)'],1);
[Grid_Wb,~,~]=plm2xyz([LM reshape(Wb',2,index/2)'],1);

% truncate the degrees out of localization range
s = find(la<Lwin+1&la>90-Lwin);

[lmcosi,~]=xyz2plm((Grid_Ht-Grid_Wt).*Gridwindow,120,'im');
Hi_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
Hi_w(s)=0;
[Grid_Hi_w,~,~]=plm2xyz([LM reshape(Hi_w',2,index/2)'],1);

[lmcosi,~]=xyz2plm(Grid_faa.*Gridwindow,120,'im');
faailm_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
faailm_w(s)=0;
[Grid_faa_w,~,~]=plm2xyz([LM reshape(faailm_w',2,index/2)'],1);

[lmcosi,~]=xyz2plm(Grid_Ht.*Gridwindow,120,'im');
Ht_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
Ht_w(s)=0;
[Grid_Ht_w,~,~]=plm2xyz([LM reshape(Ht_w',2,index/2)'],1);

[lmcosi,~]=xyz2plm(Grid_Hb.*Gridwindow,120,'im');
Hb_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
Hb_w(s)=0;
[Grid_Hb_w,~,~]=plm2xyz([LM reshape(Hb_w',2,index/2)'],1);

[lmcosi,~]=xyz2plm(Grid_Wt.*Gridwindow,120,'im');
Wt_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
Wt_w(s)=0;
[Grid_Wt_w,~,~]=plm2xyz([LM reshape(Wt_w',2,index/2)'],1);


[lmcosi,~]=xyz2plm(Grid_Wb.*Gridwindow,120,'im');
Wb_w = reshape(lmcosi(:,3:4)',1,2*length(lmcosi))';
Wb_w(s)=0;
[Grid_Wb_w,~,~]=plm2xyz([LM reshape(Wb_w',2,index/2)'],1);

% plot the gravity and flexure before localization
figure;
subplot(3,2,1)
pcolor(lonlon,latlat,Grid_Ht-Grid_Wt);shading flat;colormap jet;colorbar;xlim([200 250]);ylim([0 35])
title 'Initial surface loading'

subplot(3,2,2)
pcolor(lonlon,latlat,Grid_faa);shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Modeled FAA'

subplot(3,2,3)
pcolor(lonlon,latlat,Grid_Ht);shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Topography due to top loading Ht'

subplot(3,2,5)
pcolor(lonlon,latlat,Grid_Wt);shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Plate bending due to surface loading Wt'

subplot(3,2,4)
pcolor(lonlon,latlat,Grid_Hb);shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Topography/plate bending due to bottom loading Hb'

subplot(3,2,6)
pcolor(lonlon,latlat,Grid_Wb);shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Initial bottom loading Wb'


% plot the gravity and flexure after localization
figure;
subplot(3,2,1)
contourf(lonlon,latlat,Grid_Hi_w,'LineStyle','none');shading flat;colormap jet;colorbar;xlim([200 250]);ylim([0 35])
title 'Initial surface loading'

subplot(3,2,2)
contourf(lonlon,latlat,Grid_faa_w,'LineStyle','none');shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Modeled FAA'

subplot(3,2,3)
contourf(lonlon,latlat,Grid_Ht_w,'LineStyle','none');shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Topography due to top loading Ht'

subplot(3,2,5)
contourf(lonlon,latlat,Grid_Wt_w,'LineStyle','none');shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Plate bending due to surface loading Wt'

subplot(3,2,4)
contourf(lonlon,latlat,Grid_Hb_w,'LineStyle','none');shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Topography/plate bending due to bottom loading Hb'

subplot(3,2,6)
contourf(lonlon,latlat,Grid_Wb_w,'LineStyle','none');shading flat;colorbar;xlim([200 250]);ylim([0 35])
title 'Initial bottom loading Wb'

end



end
