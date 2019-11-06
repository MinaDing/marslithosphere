function [Tebest,Tesd,rhotbest,rhotsd,fbest,fsd] = MCMC_flexure_VolcanicMontes_coarsened(gamma,Te_range,rhot_range,f_range,Te_isd,rhot_isd,f_isd,Te_ini,rhot_ini,f_ini,N,burnin,DEGS,FADM,FADM_sd,capid,degmin,degmax)
% input: 
%   gamma: coarsening parameter
%   Te_range, rhot_range, f_range: A priori range for the parameters
%   Te_isd, rhot_isd, ...: initial guess of the standard deviation 
%   Te_ini, rhot_ini,...: initial values for the Markov chain
%   N: length of the Markov chain 
%   burnin: burnin period (may adjust if not ideal in the first run)
%   DEGS,FADM,FADM_sd: input spectral data to fit
%   capid: the capid necessary for forward modeling
%   degmin,degmax: min and max degree for inversion


ifplot = 0;
format long
%disp 'Apply MCMC method to constrain flexural parameters, Te, rhot, and f'  
degs_ix = find(DEGS>=degmin&DEGS<=degmax);
%disp 'Fit free-air admittance'

% Bayesian parameter estimation
acc=0; % number of accepted step after burnin steps
dim=3; % total # of parameters: Te, rhot, f
prop_Var=[Te_isd rhot_isd f_isd]; % pre-assumed Gaussian variance

% Initialize Markov chain
X = zeros(dim,N); 
X(:,1) = [Te_ini;rhot_ini;f_ini]; % initial values

%disp('Initial setup for parameter steps - SD(Te), SD(rhot), SD(f)')
%disp(num2str(prop_Var))

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
f = X(3,1);

[fadm,~] = mcgovern_synthetic_admcor(Te,rhot,f,capid,DEGS);

    if ifplot
    h1 = figure;
    errorbar(DEGS,FADM,FADM_sd,'k');hold on;
%    plot(DEGS,FADM,'k');hold on;    
    plot(DEGS,fadm,'r')
    ylabel('Free-air adm. (mGal/km)')
    set(gca,'Ytick',-100:50:260,'ylim',[-100 300],'Xtick',0:20:120,'xlim',[min(DEGS) max(DEGS)],'YMinorTick','on','XMinorTick','on')
    title(['k=1,accepted steps=0,Te=' num2str(X(1,1)) ' km, rhot=' num2str(X(2,1)) ' kg/m^3, f=' num2str(X(3,1))]);
    disp(['k=1,accepted steps=0,Te=' num2str(X(1,1)) ' km, rhot=' num2str(X(2,1)) ' kg/m^3, f=' num2str(X(3,1))]);
    
    h2 = figure;
    
    end

% Calculate prior 
currPrior = unifpdf(X(1,1),Te_range(1),Te_range(2))*unifpdf(X(2,1),rhot_range(1),rhot_range(2))*unifpdf(X(3,1),f_range(1),f_range(2)); % continuous uniform pdf
%currPrior = unifpdf(X(1,1),Te_range(1),Te_range(2))*unifpdf(X(2,1),rhot_range(1),rhot_range(2))*unifpdf(X(3,1),f_range(1),f_range(2))*normpdf(X(4,1),0,0.1); % continuous uniform pdf

% Calculate likelihood
%currLikely = prod(normpdf(FADM,fadm,sqrt(Sigma2e)))*prod(normpdf(FCOR,fcor,sqrt(Sigma2e1)))*prod(normpdf(BCOR,bcor,sqrt(Sigma2e1))); % multi-variate Gaussian;log current likelihood
%LcurrLikely = log(prod(normpdf(FADM,fadm,FADM_sd)))+log(prod(normpdf(FCOR,fcor,FCOR_sd)))+log(prod(normpdf(BCOR,bcor,BCOR_sd))); % multi-variate Gaussian;log current likelihood
%LcurrLikely = log(prod(normpdf(FADM,fadm,FADM_sd)))+log(prod(normpdf(FCOR,fcor,FCOR_sd))); % multi-variate Gaussian;log current likelihood
LcurrLikely = log(prod(normpdf(FADM(degs_ix),fadm(degs_ix),FADM_sd(degs_ix)))); 
LcurrLikely = 1/gamma*LcurrLikely;


% Initialize the results
fadm_best = fadm*0;
fadm_sd = fadm_best;

% Construct Markov chain    
for k=2:N
    Yprop=X(:,k-1) + (randn(1,dim)*R)'; % random walking   
    % calculate the synthetic curves 
    Te = Yprop(1)*1e3;
    rhot = Yprop(2);        
    f = Yprop(3);
    [fadm,~]=mcgovern_synthetic_admcor(Te,rhot,f,capid,DEGS);
    
    % Update prior and likelihood
    propPrior = unifpdf(Yprop(1),Te_range(1),Te_range(2))*unifpdf(Yprop(2),rhot_range(1),rhot_range(2))*unifpdf(Yprop(3),f_range(1),f_range(2)); % continuous uniform pdf
    LpropLikely = log(prod(normpdf(FADM(degs_ix),fadm(degs_ix),FADM_sd(degs_ix))));
    LpropLikely = 1/gamma*LpropLikely;          
    
    % Calculate acceptance probability
    if (currPrior==0||LcurrLikely==-Inf)
        X(:,k) = X(:,k-1);
    else
        a = min(1,exp(LpropLikely)/exp(LcurrLikely)*propPrior/currPrior);
        if (rand()<=a)
            X(:,k) = Yprop;
            LcurrLikely = LpropLikely;
            currPrior = propPrior;
            fadm_curr = fadm;
            %disp 'Markov chain updated!'
    %%
    if ifplot
    figure(h1);
    cla
    errorbar(DEGS,FADM,FADM_sd,'k');hold on;
%    plot(DEGS,FADM,'k');hold on;    
    plot(DEGS,fadm,'r')
    ylabel('Free-air adm. (mGal/km)')
    set(gca,'Ytick',-100:50:260,'ylim',[-100 300],'Xtick',0:20:120,'xlim',[min(DEGS) max(DEGS)],'YMinorTick','on','XMinorTick','on')
    title(['k=' num2str(k) ',accepted steps=' num2str(acc) ',Te=' num2str(X(1,k)) ' km, rhot=' num2str(X(2,k)) ' kg/m^3, f=' num2str(X(3,k))]);     
    pause(0.001)
    end
    
        else
            X(:,k) = X(:,k-1);
        end
    end
    
    if k==burnin
        fadm_min = fadm_curr;
        fadm_max = fadm_curr;
    end
    if k>burnin
        fadm_min = min(fadm_min,fadm_curr);
        fadm_max = max(fadm_max,fadm_curr);
    end
    if(k>=burnin)
        fadm_best = fadm_best+fadm_curr;
        fadm_sd = fadm_sd+fadm_curr.^2;
        acc = acc+1;
    end
    
 %   disp(['k=' num2str(k) ',accepted steps=' num2str(acc) ',Te=' num2str(X(1,k)) ' km, rhot=' num2str(X(2,k)) ' kg/m^3, f=' num2str(X(3,k))]);

    % Update the adaptive variance matrix
    if mod(k,k0)==1
        V(:,:,k)=sp*cov(X(:,1:k-1)')+epsilon*eye(dim,dim);
        R=chol(V(:,:,k));
    else
        V(:,:,k)=V(:,:,k-1);
    end
    
    if ifplot
    figure(h2);
    subplot1(3,1);
    subplot1(1)
    cla
    plot(X(1,:),'k');hold on;
    line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
    ylabel 'Te (km)'
    
    subplot1(2)
    cla
    plot(X(2,:),'k');hold on;
    line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
    ylabel '\rhot (kg/m^3)'

    subplot1(3)
    cla
    plot(X(3,:),'k');hold on;
    line([burnin burnin],get(gca,'YLim'),'color','k','linestyle','--')
    ylabel 'f'
    pause(0.001)
    end
    
  %  if mod(k,100)==0 % save the results after generating every 100 samples
  %      save tmp.mat
  %  end
%     
    
end

fadm_best = fadm_best/acc;
fadm_sd = sqrt(fadm_sd/acc-fadm_best.^2);
Tebest = mean(X(1,burnin:N));
Tesd = std(X(1,burnin:N));
fbest = mean(X(3,burnin:N));
fsd = std(X(3,burnin:N));
rhotbest = mean(X(2,burnin:N));
rhotsd = std(X(2,burnin:N));

%save tmp.mat
save(['Volcano10_Olympus_3para_' num2str(degmin) num2str(degmax) '_' num2str(gamma) '.mat']);

end