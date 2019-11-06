function MCMC_ImpactSpatial_coarsened(gamma,Te_range,rhot_range,Te_isd,rhot_isd,Te_ini,rhot_ini,N,burnin,Dvec,FAAvec,FAAvec_sd,la,LM,topox_lm,Grid_d,TOPO_end)

ifplot = 0;
% Apply MCMC method to lithospheric flexure to solve for three parameters
% including: Te (lithospheric thickness); rhot (surface load density); f
% (sub-surface/surface load ratio)

% Input: 
% Te_range: prior range for Te (km)
% rhot_range: prior range for rhot
% f_range: prior range for f

% Te_isd: Intial standard deviation for Te to be upated using adaptive vairance algorithm
% rhot_isd: Initial standard deviation for rhot
% f_isd: Initial s.d. for f
% alpha_isd: Initial s.d. for alpha

% N: assumed number of iterations
% burnin: burnin period to adjust manually

% flexure_admcor2: a function to calculate synthetic admittance correlation curves, constant parameters in sh_admcor_ini.mat
% e.g. [fadm,fcor,badm,bcor]=flexure_admcor2(Te,f*rhot/rhob,alpha,DEGS,rhot,rhob,rhoc,rhom,E,nu,Rplanet,gzero,Hc,Robs);

% Written by Chen Gu (guchch@mit.edu) Jan 4, 2016
% Revised by Min Ding (dingmin@mit.edu) April 5, 2016

format long
%addpath('/Users/dingmin/Documents/MATLAB/Mars','-end');
%addpath('/Users/dingmin/Documents/MATLAB/Marsflexure','-end');

disp 'Apply MCMC method to constraining flexural parameters, Te, rhot, and f'  

% number of iterations
if ifplot
figure;
subplot1(3,1)
end

%% Fit all the three curves

% Bayesian parameter estimation
acc=0; % number of accepted step after burnin steps

dim=2; % total # of parameters: Te, rhot, f
prop_Var=[Te_isd rhot_isd]; % pre-assumed Gaussian variance

% Initialize Markov chain
X = zeros(dim,N); 
X(:,1) = [Te_ini;rhot_ini]; % initial values

disp('Initial setup for parameter steps - SD(Te), SD(rhot), SD(f)')
disp(num2str(prop_Var))

% Adaptive variance
sp=2.38^2/dim; % constant 
epsilon=0.1; % 0-0.1 Matrix diagnal
V = zeros(dim,dim,N); % covariance matrix
V(:,:,1) = diag(prop_Var); % initial
R = chol(V(:,:,1)); % decompose ... Lower/upper matrix
k0 = 50; % after k0 steps, covariance matrix is updated; k0 can be varied if needed
%R = diag(prop_Var);

% Calculate synethic curves 
Te = X(1,1)*1e3;
rhot = X(2,1);        

% the localization depends on the beta value 
faa = GravityPredict(rhot,Te,Dvec,la,LM,topox_lm,Grid_d,TOPO_end);
faa = faa-faa(1)+FAAvec(1);

% Calculate prior 
currPrior = unifpdf(X(1,1),Te_range(1),Te_range(2))*unifpdf(X(2,1),rhot_range(1),rhot_range(2)); % continuous uniform pdf

% Calculate likelihood
%currLikely = prod(normpdf(FADM,fadm,sqrt(Sigma2e)))*prod(normpdf(FCOR,fcor,sqrt(Sigma2e1)))*prod(normpdf(BCOR,bcor,sqrt(Sigma2e1))); % multi-variate Gaussian;log current likelihood
LcurrLikely = log(prod(normpdf(FAAvec',faa,FAAvec_sd))); % multi-variate Gaussian;log current likelihood
LcurrLikely = 1/gamma*LcurrLikely;


% Initialize the results
faa_best = faa*0;
faa_sd = faa_best;

% Construct Markov chain    
for k=2:N

    Yprop=X(:,k-1) + (randn(1,dim)*R)'; % random walking
        
    % calculate the synthetic curves 
    Te = Yprop(1)*1e3;
    rhot = Yprop(2);  
    faa = GravityPredict(rhot,Te,Dvec,la,LM,topox_lm,Grid_d,TOPO_end);
    faa = faa-faa(1)+FAAvec(1);
    
    % Update prior and likelihood
    propPrior = unifpdf(Yprop(1),Te_range(1),Te_range(2))*unifpdf(Yprop(2),rhot_range(1),rhot_range(2)); % continuous uniform pdf
    LpropLikely = log(prod(normpdf(FAAvec',faa,FAAvec_sd)));% multi-variate Gaussian of Yi^{k}
    LpropLikely = 1/gamma*LpropLikely;
    
    % Calculate acceptance probability
    %if (currPrior==0||currLikely==0)
    if (currPrior==0||LcurrLikely==-Inf)
        X(:,k) = X(:,k-1);
    else
        %a = min(1,propLikely/currLikely*propPrior/currPrior);
        a = min(1,exp(LpropLikely)/exp(LcurrLikely)*propPrior/currPrior);
        if (rand()<=a)
            X(:,k) = Yprop;
            LcurrLikely = LpropLikely;
            currPrior = propPrior;
            faa_curr = faa;
            if ifplot
                subplot1(1);
                cla
                errorbar(DEGS,FADM,FADM_sd,'k');hold on;
            %    plot(DEGS,FADM,'k');hold on;    
                plot(DEGS,fadm,'r')
                ylabel('Free-air adm. (mGal/km)')
                set(gca,'Ytick',-100:50:260,'ylim',[-100 300],'Xtick',0:20:120,'xlim',[min(DEGS) max(DEGS)],'YMinorTick','on','XMinorTick','on')
                subplot1(2);
                cla
                errorbar(DEGS,FCOR,FCOR_sd,'k');hold on;
            %    plot(DEGS,FCOR,'k');hold on;
                plot(DEGS,real(fcor),'r')
                ylabel('Free-air cor.')
                set(gca,'Ytick',-0.5:0.5:1,'Ylim',[-1 1],'Xtick',0:20:120,'xlim',[min(DEGS) max(DEGS)],'YMinorTick','on','XMinorTick','on')


                subplot1(3);
                cla
                errorbar(DEGS,BCOR,BCOR_sd,'k');hold on;
            %    plot(DEGS,BCOR,'k');hold on;
                plot(DEGS,real(bcor),'r')
                ylabel('Bouguer cor.')
                set(gca,'Ytick',-1:0.5:1,'Ylim',[-1 1],'Xtick',0:20:120,'xlim',[min(DEGS) max(DEGS)],'YMinorTick','on','XMinorTick','on')
                pause(0.01)      
            end
    
        else
            X(:,k) = X(:,k-1);
        end
    end
    
    if(k>=burnin)
        faa_best = faa_best+faa_curr;
        faa_sd = faa_sd+faa_curr.^2;
        acc = acc+1;
    end
    
    
    % Update the adaptive variance matrix
    if mod(k,k0)==1
        V(:,:,k)=sp*cov(X(:,1:k-1)')+epsilon*eye(dim,dim);
        R=chol(V(:,:,k));
    else
        V(:,:,k)=V(:,:,k-1);
    end
    
    disp(['k=' num2str(k) ',accepted steps=' num2str(acc) ',Te=' num2str(X(1,k)) ' km, rhot=' num2str(X(2,k)) ' kg/m^3']);


end

faa_best = faa_best/acc;
faa_sd = sqrt(faa_sd/acc-faa_best.^2);


Tebest = mean(X(1,burnin:N));
Tesd = std(X(1,burnin:N));
rhotbest = mean(X(2,burnin:N));
rhotsd = std(X(2,burnin:N));

save tmp.mat

end