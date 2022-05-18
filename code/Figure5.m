% Bayesian model simulations of sound discrimination for Read Lab
% Sashank Pisupati 02/20/20


%% Model specifications
clear all
plotSim = 1; % Plot simulations from ideal observer
fitting = 1; % Fit rodent data
plotFit = 1; % Plot data fits

% -----------------------Fitting specifications---------------------------

% Models to be compared
priorTypes = {'onset'};     % Empirical priors constructed from onset/offset/both slopes of vocalization data
lapseTypes = {'fixed','inattention'};       % Models of lapses - fixed (random guesses) or inattentive (guesses from prior)
noiseTypes = {'fixed','variable'};          % Duration noise constraint - fixed or varying across slope conditions
integrationTypes = {'perfect','imperfect'}; % Integration constraint - perfect integration for hold duration, or imperfect (tInt<tHold)

% Load data, empirical parameters
try
%     load('dataset_Smith_Jane.mat')% Across-subject across-session pooled data for duration discrimination task
    load('/data/dataset_Jane_Pisupati.mat')
catch
    fprintf('\nCannot find data directory..proceeding with simulated data\n')
end
if exist('dataset','var')
    datacat = [dataset(:).data];
    nTrials = sum([datacat.n]);                 % Total number of trials in dataset
    tHold = [dataset.relativeHoldTime];         % Empirical median hold times in units of bursts ~ integration times
else
    nTrials = 30000;
    tHold = [1 1.4 2.1 1 1.4 2.1];
end

% -------------------------Simulation specifications----------------------

priorType = 'onset';
lapseType = 'fixed';
noiseType = 'variable';
integrationType = 'perfect';

% Generative parameters
switch noiseType
    case 'fixed'
        noiseSim = [50,50,50,50,50,50];% Generative sensory noise parameter
    case 'variable'
        noiseSim = [44,44,44,30,30,30];% Generative sensory noise parameter
end
switch integrationType
    case 'perfect'
        tInt = tHold;                  % Perfect integration for hold duration
    case 'imperfect'
        tInt = 0.9*tHold;              % Imperfect integration times (tInt<tHold)
end
sigmaSim = noiseSim./sqrt(tHold);      % Uncertainty reduction through integration
lambdaSim = 0.25;                       % Generative lapse parameter


% --------------------------Plotting specifications-----------------------

tasks = {'plateauSlowOnOff','plateauSlowOnOff','plateauSlowOnOff','plateauFastOnOff','plateauFastOnOff','plateauFastOnOff'};
holdTimes = {'200ms single','200ms multi','600ms','200ms single','200ms multi','600ms'};
colors = flip({[0,0,142] ./ 255, [7,72,255] ./ 255, [95,177,255] ./ 255, [182,0,0] ./ 255,[255,7,45] ./ 255,[255,95,146] ./ 255});
priorInds = [2,2,2,2,2,2];
perfInds = [4,5,6,4,5,6];
markerSize = 20;
lineWidth = 2;

%% Simulate synthetic data based on empirical parameters (Figure 2)
if plotSim
    
    % Plot empirical prior
    jointPrior = getJointPrior(priorType);
    figure(1);clf;set(gcf,'color','w')
    subplot(2,3,1)
    x = [0:50:2000];
    y = [0:50:800];
%     y = [0: 0.001: 0.6];
    [X,Y] =meshgrid(x,y);
    z = log(mvnpdf([reshape(X,numel(X),1),reshape((Y),numel(Y),1)],jointPrior.Mu,jointPrior.Sigma));
    contour(x,y,reshape(z,size(X,1),size(X,2)),15)
    colormap('winter')
    
    % Plot simulated data
    for i = 1:length(tasks)
        
        % Get appropriate stimulus, prior & parameters for current task
        task = tasks{i};
        [s,prior] = getStim(task,jointPrior);
        stim = s.(s.taskDim);      % Stimulus set
        params.muP = prior.mu;     % Prior mean
        params.sigP = prior.sig;   % Prior standard deviation
        rules.mu0 = s.mu0;         % Category boundary
        params.sigS = sigmaSim(i); % Sensory noise parameter
        params.lambda = lambdaSim; % Lapse rate parameter
        
        % Plot stimulus design for current task
        subplot(2,3,1)
        hold on
        plot(s.plateauDur,s.onOffFc,'o','MarkerFaceColor',colors{i},'color',colors{i},'MarkerSize',5)
        ylabel('Onset-Offset Slope (A/s)')
        xlabel('Plateau duration (ms)')
        ylim([y(1),y(end)])
        xlim([x(1),x(end)])
        title('Task design')
        
        % Plot duration priors for current task
        subplot(2,3,priorInds(i))
        xPrior = [0:(1500)/99:1500];
        yMaxPrior = 0.005;
        xlabel('Plateau duration (ms)')
        title([ priorType, ' prior'])
        plot([prior.mu,prior.mu],[0,yMaxPrior],'--','color',colors{i})
        plot([s.mu0,s.mu0],[0,yMaxPrior],'k-')
        hold on
        plot(xPrior,normpdf(xPrior,prior.mu,prior.sig),'color',colors{i},'LineWidth',lineWidth)
        ylim([0,yMaxPrior]);
        xlim([xPrior(1),xPrior(end)])
        
        % Simulate synthetic data
        thisStim = randsample(stim,round(nTrials/length(tasks)),true);
        thisChoice = simulateIdealObserver(params,rules,thisStim,lapseType);
        allStim{i} = thisStim;
        allChoices{i} = thisChoice;
        
        % Generate & plot psychometric curves
        subplot(2,3,perfInds(i));
        hold on
        clear nHighResp nResp pChoice
        for j  = 1:length(stim)
            % Get fraction chose "High" response
            nHighResp(j) = sum(thisChoice(thisStim==stim(j)));
            nResp(j) = sum(thisStim==stim(j));
            pChoice(j) = nHighResp(j)/nResp(j);
            % Get 95% Wilson binomial confidence intervals
            z = 1.96;
            ciUpper(j) = (pChoice(j) + z^2/(2*nResp(j)) + z * sqrt(pChoice(j)*(1-pChoice(j))/nResp(j) + z^2/(4*nResp(j)^2))) / (1 + z^2/nResp(j));
            ciLower(j) = (pChoice(j) + z^2/(2*nResp(j)) - z * sqrt(pChoice(j)*(1-pChoice(j))/nResp(j) + z^2/(4*nResp(j)^2))) / (1 + z^2/nResp(j));
            % Plot confidence intervals
%             plot([stim(j),stim(j)],[ciUpper(j),ciLower(j)],'color',colors{i});
        end
        % Plot datapoints
%         scatter(stim, pChoice,'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i},'sizeData',markerSize);
        
        %Plot predicted psychometric fit
        xFit = [0.8*stim(1):(1.1*stim(end)-0.8*stim(1))/99:1.1*stim(end)];
        yFit = pChoiceIdeal([params.sigS,params.muP,params.sigP,params.lambda],rules,xFit,lapseType);
        plot(xFit,yFit,'color',colors{i},'LineWidth',lineWidth)
        plot([rules.mu0,rules.mu0],[0,1],'k-')
        xlim([xFit(1),xFit(end)])
        ylim([0,1])
        title(['Simulation - ',holdTimes{i}]);
        
        % Plot predicted bias parameter
        wP = 1/prior.sig^2;
        wS = 1/params.sigS^2;
        wTot = wP+wS;
        biasSim(i) = rules.mu0 -(wS*rules.mu0 + wP*prior.mu)/wTot;
        plot(rules.mu0+[biasSim(i),biasSim(i)],[0,1],'--','color',colors{i})
    end
end

saveas(gcf,'/results/model_simulation.pdf')

%% Joint fits to rodent data (if available - otherwise uses simulated data from above)

% -----------------Get data----------------
if exist('dataset','var')
    for i = 1:length(tasks)
        % Summarized data
        clear stim nHighResp nResp pChoice
        validResp = dataset(i).data.y~=0;
        stim = 100*dataset(i).data.x(validResp);
        nHighResp = dataset(i).data.y(validResp);
        nResp = dataset(i).data.n(validResp);
        
        % Expand into (unordered) trials
        thisStim = zeros(1,sum(nResp));
        thisChoice = zeros(1,sum(nResp));
        j0 = 0;
        for j = 1:length(stim)
            thisStim(j0+1:j0+nResp(j))=stim(j);
            thisChoice(j0+1:j0+nHighResp(j))=1;
            j0 = j0+nResp(j);
        end
        
        allStim{i} = thisStim;
        allChoices{i} = thisChoice;
    end
end

%----------------Fit models-------------------

if fitting
    k = 0;
    % Prior specification
    for p = 1:length(priorTypes)
        
        for i = 1:length(tasks)
            task = tasks{i};
            priorType = priorTypes{p};
            jointPrior = getJointPrior(priorType); % Empirical joint prior of durations & slopes
            [s,prior] = getStim(task,jointPrior);
            rules.mu0 = s.mu0;
            allPriors{i} = prior;
            allRules{i} = rules;
        end
        
        
        % Sensory noise specification
        for q = 1:length(noiseTypes)
            noiseType = noiseTypes{q};
            % Integration model
            for r = 1:length(integrationTypes)
                integrationType = integrationTypes{r};
                % Lapse model
                for s = 1:length(lapseTypes)
                    lapseType = lapseTypes{s};
                    k = k+1;
                    fprintf('Fitting model %d/8\n',k)
                    clear thetaJoint thetaJointEst
                    
                    % ------------------Define constraints---------------------
                    % Same lapse rate across conditions
                    Aeq =  [0,1,0,-1,0,0,0,0,0,0,0,0;...
                        0,0,0,1,0,-1,0,0,0,0,0,0;...
                        0,0,0,0,0,1,0,-1,0,0,0,0;...
                        0,0,0,0,0,0,0,1,0,-1,0,0;...
                        0,0,0,0,0,0,0,0,0,1,0,-1];
                    Beq = [0;0;0;0;0];
                    
                    switch integrationType
                        case 'perfect'
                            % Sensory noise goes down with sqrt of integration time
                            Aeq = [Aeq;1,0,-sqrt(tHold(2)),0,0,0,0,0,0,0,0,0;...
                                1,0,0,0,-sqrt(tHold(3)),0,0,0,0,0,0,0;...
                                0,0,0,0,0,0,1,0,-sqrt(tHold(5)),0,0,0;...
                                0,0,0,0,0,0,1,0,0,0,-sqrt(tHold(6)),0];
                            Beq = [Beq;0;0;0;0];
                            nlcon = [];
                            switch noiseType
                                case 'fixed'
                                    % Same sensory noise across slope conditions
                                    Aeq = [Aeq;1,0,0,0,0,0,-1,0,0,0,0,0];
                                    Beq = [Beq;0];
                                    nParams = 2;
                                case 'variable'
                                    % Each slope condition get own sensory noise
                                    nParams = 3;
                            end
                        case 'imperfect'
                            % Integration times t1 & t2 are free parameters
                            switch noiseType
                                case 'fixed'
                                    % Same sensory noise across slope conditions
                                    Aeq = [Aeq;1,0,0,0,0,0,-1,0,0,0,0,0;...
                                        0,0,1,0,0,0,0,0,-1,0,0,0;...
                                        0,0,0,0,1,0,0,0,0,0,-1,0];
                                    Beq = [Beq;0;0;0];
                                    nlcon = [];
                                    nParams = 4;
                                case 'variable'
                                    % Variabl sensory noise subject to integration
                                    nlcon = @(x)integration_constraint(x);
                                    nParams = 5;
                            end
                    end
                    
                    % ---------------------Fit parameters-------------------
                    opts = optimoptions('fmincon');
                    opts.Display = 'none';
                    LB = [0,0,0,0,0,0,0,0,0,0,0,0];
                    UB = [Inf,1,Inf,1,Inf,1,Inf,1,Inf,1,Inf,1];
                    paramRange = [100,1,100,1,100,1,100,1,100,1,100,1];
                    theta0 = LB + paramRange.*rand(1,length(LB));
                    [thetaEst,nll] = ...
                        fmincon(@(theta)nllJoint(theta,...
                        allPriors,allRules,allStim,allChoices,lapseType),...
                        theta0,[],[],Aeq,Beq,LB,UB,nlcon,opts);
                    

                    
                    t1 = (thetaEst(1)/thetaEst(3))^2;
                    t2 = (thetaEst(1)/thetaEst(5))^2;
                    
                    fits(k).priorType = priorType;
                    fits(k).noiseType = noiseType;
                    fits(k).integrationType = integrationType;
                    fits(k).lapseType = lapseType;
                    fits(k).nll = nll;
                    fits(k).bic = nParams*log(nTrials) + 2*nll;
                    fits(k).aic = 2*nParams + 2*nll;
                    fits(k).theta = thetaEst;
                    fits(k).t1 = t1;
                    fits(k).t2 = t2;
                end
            end
        end
    end
else
    try
        load('factorial_model_fits.mat')
    catch
        fprintf('\nCannot find directory with model fits! Please turn on fitting\n')
    end
end

[~,ind_best_BIC] = min([fits.bic]-min([fits.bic]));
[~,ind_best_AIC] = min([fits.aic]-min([fits.aic]));

%% Plot fits

%----------------Psychometric curve fits--------------------
if plotFit
    
    % Get model specification
    ind = ind_best_BIC;
    %  ind = ind_best_AIC;
    priorType = fits(ind).priorType;
    noiseType = fits(ind).noiseType;
    integrationType = fits(ind).integrationType;
    lapseType = fits(ind).lapseType;
    thetaEst = fits(ind).theta;
    [jointPrior] = getJointPrior(priorType);
    
    hfig = figure(2);clf;set(gcf,'color','w')
    
    % Plot empirical prior
    subplot(2,3,1)
%     x = [0:50:400];
%     y = [0:200:2000];
    [X,Y] =meshgrid(x,y);
    z = log(mvnpdf([reshape(X,numel(X),1),reshape((Y),numel(Y),1)],jointPrior.Mu,jointPrior.Sigma));
    colormap('winter')
    contour(x,y,reshape(z,size(X,1),size(X,2)))
    
    % Plot data split by condition
    for i = 1:length(tasks)
        task = tasks{i};
        [s,prior] = getStim(task,jointPrior);
        rules.mu0 = s.mu0;
        
        % Plot stimulus design
        subplot(2,3,1)
        hold on
        plot(s.plateauDur,s.onOffFc,'o','MarkerFaceColor',colors{i},'color',colors{i},'MarkerSize',5)
        ylabel('Onset-Offset Slope (log dB/s)')
        xlabel('Plateau duration (ms)')
        ylim([y(1),y(end)])
        xlim([x(1),x(end)])
        title('Task design')
        
        % Plot conditional priors
        subplot(2,3,priorInds(i))
        hold on
        xFit = [0:(1500)/99:1500];
        xlabel('Plateau duration (ms)')
        yMaxPrior = 0.005;
        title([ priorType, ' prior'])
        ylim([0,yMaxPrior]);
        plot([prior.mu,prior.mu],[0,yMaxPrior],'--','color',colors{i})
        plot([s.mu0,s.mu0],[0,yMaxPrior],'k-')
        plot(xFit,normpdf(xFit,prior.mu,prior.sig),'color',colors{i},'LineWidth',lineWidth)
        xlim([xFit(1),xFit(end)])
        
        
        % Generate psychometric curves
        subplot(2,3,perfInds(i));
        hold on
        title(['Rat data - ',holdTimes{i}])
        thisStim = allStim{i};
        thisChoice = allChoices{i};
        
        clear nHighResp nResp pChoice
        for j  = 1:length(stim)
            % Get fraction chose "High" response
            stim = unique(thisStim);
            nHighResp(j) = sum(thisChoice(thisStim==stim(j)));
            nResp(j) = sum(thisStim==stim(j));
            pChoice(j) = nHighResp(j)/nResp(j);
            % Get 95% Wilson binomial confidence intervals
            z = 1.96;
            ciUpper(j) = (pChoice(j) + z^2/(2*nResp(j)) + z * sqrt(pChoice(j)*(1-pChoice(j))/nResp(j) + z^2/(4*nResp(j)^2))) / (1 + z^2/nResp(j));
            ciLower(j) = (pChoice(j) + z^2/(2*nResp(j)) - z * sqrt(pChoice(j)*(1-pChoice(j))/nResp(j) + z^2/(4*nResp(j)^2))) / (1 + z^2/nResp(j));
            % Plot confidence intervals
            plot([stim(j),stim(j)],[ciUpper(j),ciLower(j)],'color',colors{i});
        end
        % Plot datapoints
        scatter(stim, pChoice,'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i},'sizeData',markerSize);
        
        % Plot 95% Wilson binomial confidence intervals
        z = 1.96;
        ciUpper = (pChoice + z^2./(2*nResp) + z .* sqrt(pChoice.*(1-pChoice)./nResp + z^2./(4*nResp.^2))) ./ (1 + z^2./nResp);
        ciLower = (pChoice + z^2./(2*nResp) - z .* sqrt(pChoice.*(1-pChoice)./nResp + z^2./(4*nResp.^2))) ./ (1 + z^2./nResp);
        for j  = 1:length(stim)
            plot([stim(j),stim(j)],[ciUpper(j),ciLower(j)],'color',colors{i});
        end
        plot([rules.mu0,rules.mu0],[0,1],'k-')
        
        
        % Plot fits
        sigmaFit = thetaEst(2*i-1);
        lambdaFit = thetaEst(2*i);
        xFit = [0.8*stim(1):(1.1*stim(end)-0.8*stim(1))/99:1.1*stim(end)];
        yFit = pChoiceIdeal([sigmaFit,prior.mu,prior.sig,lambdaFit],rules,xFit,lapseType);
        plot(xFit,yFit,'-','LineWidth',1,'color',colors{i})
        xlim([xFit(1),xFit(end)])
        ylim([0,1])
        
        % Extract fit bias parameter, bootstrap SE
        wP = 1/prior.sig^2;
        wS = 1/sigmaFit^2;
        wTot = wP+wS;
        biasFit(i) = (wS*rules.mu0 + wP*prior.mu)/wTot-rules.mu0;
        
    end
    
    
% ------------Plot data vs fit biases as a function of hold times----------
    if exist('dataset','var')
        for i = 1:length(dataset)
            biasData(i) = rules.mu0-dataset(i).params_palamedes(1);
            biasSEData(i) = dataset(i).paramsSE_palamedes(1);
        end
    else
        biasData = -biasSim;
        biasSEData = zeros(1,6);
    end
    
    t1 = fits(ind).t1;
    t2 = fits(ind).t2;
    tInt = [1 t1 t2];
    
    subplot(2,3,3)
    hold on
    h(1) = plot(tInt,biasData(1:3)-biasData(4:6),'kx--','LineWidth',1);
    h(2) = plot(tInt,biasFit(1:3)-biasFit(4:6),'bo-','LineWidth',1);
    for t = 1:3
        dataSE = sqrt(biasSEData(t)^2 + biasSEData(t+3)^2);
        plot([tInt(t),tInt(t)],[biasData(t)-biasData(t+3)+dataSE, biasData(t)-biasData(t+3)-dataSE],'k-')
    end
    ylim([0,40])
    xlim([0.5,2.5])
    box off
    ylabel('\Delta Bias (Slow-Fast)')
    xlabel('Median Hold Time (s)')
    legend(h,'Data','Fit')
    
end

saveas(gcf,'/results/figure_5_model_fits.pdf')

%----------------------Plot factorial model comparison----------------------
figure(3);clf;set(gcf,'color','w')
    % Prior specification
    ind = 0;
    for p = 1:length(priorTypes)
        priorType = priorTypes{p};
        k = 0;
        % Sensory noise specification
        for q = 1:length(noiseTypes)
            noiseType = noiseTypes{q};
            % Integration model
            for r = 1:length(integrationTypes)
                integrationType = integrationTypes{r};
                % Lapse model
                for s = 1:length(lapseTypes)
                    lapseType = lapseTypes{s};
                    ind = ind+1;
                    k = k+1;
                    tag{k} = [noiseType(1:3),'-',integrationType(1:4),'-',lapseType(1:3)];
                    bic(k) = fits(ind).bic;
                    aic(k) = fits(ind).aic;
                end
            end
        end
    end
    bic = (bic-min(min(bic)));
    aic = (aic-min(min(aic)));
    subplot(1,2,1)
    bar(log(bic))
%     legend({'onset','offset','both'},'box','off')
    set(gca,'XTick',[1:8],'XTickLabel',tag,'XTickLabelRotation',45)
    axis square
    box off
    ylabel(' log \Delta BIC')
    title('BIC')
    subplot(1,2,2)
    bar(log(aic))
%     legend({'onset','offset','both'},'box','off')
    set(gca,'XTick',[1:8],'XTickLabel',tag,'XTickLabelRotation',45)
    axis square
    box off
    ylabel(' log \Delta AIC')
    title('AIC')
    set(gcf,'Units','inches');

saveas(gcf,'/results/figure_5_model_comparison.pdf')

    
%%
% figure(4)%
% zmax = 5000;
% subplot(2,2,1)
% imagesc(([bic([1,3]);bic([5,7])]))
% caxis(log([1,zmax]))
% colorbar
% subplot(2,2,2)
% imagesc(([bic([2,4]);bic([6,8])]))
% caxis(log([1,zmax]))
% colorbar
% subplot(2,2,3)
% imagesc(([aic([1,3]);aic([5,7])]))
% caxis(log([1,zmax]))
% colorbar
% subplot(2,2,4)
% imagesc(([aic([2,4]);aic([6,8])]))
% colormap('pink')
% caxis(log([1,zmax]))
% colorbar
%% Functions
%--------------------------Get stimulus set used in task------------------
function [stim,condPrior] = getStim(task,jointPrior)
switch task
    case 'plateauSlowOnOff'
        stim.plateauDur = [100,130,160,175,190,220,250]; %ms
        stim.mu0 = 175;
%         stim.onOffFc = 255*ones(1,7); % dB/s
        stim.onOffFc = 83.7216437339341*ones(1,7); % dB/s
%         stim.onOffFc = .01465129 * ones (1,7);
        stim.taskDim = 'plateauDur';
        [condPrior.mu,condPrior.sig] = getCondPrior(jointPrior.Mu, jointPrior.Sigma, 2, stim.onOffFc(1));
        
    case 'plateauFastOnOff'
        stim.plateauDur = [100,130,160,175,190,220,250]; %ms
        stim.mu0 = 175;
%         stim.onOffFc = 1636*ones(1,7); % dB/s
        stim.onOffFc = 535.816943983612*ones(1,7); % dB/s
%         stim.onOffFc = .09376797 * ones (1,7);
        stim.taskDim = 'plateauDur';
        [condPrior.mu,condPrior.sig] = getCondPrior(jointPrior.Mu, jointPrior.Sigma, 2, stim.onOffFc(1));
end
end

%---------------Ideal observer with no/fixed/inattentive lapses------------
function pChoice = pChoiceIdeal(params, rules, stim, lapseType)
% Sensory noise
sigS = params(1);
% Prior over stimuli
muP = params(2);
sigP = params(3);
% Lapse rate
lambda = params(4);
lambda0 = 1e-4; % Small additional baseline lapse rate of 1e-4 for stability
% True criterion (assumed to be known)
mu0 = rules.mu0;

% Weights on prior, likelihood
wP = 1/sigP^2;
wS = 1/sigS^2;
wTot = wP+wS;
% Mean, SD of posterior probability of stimulus given noisy observation
mu = (wS*stim + wP*muP)/wTot;
sig = sqrt(wS)/wTot;
% Posterior probability of a "higher than" choice, marginalized over noisy
% observations of stimulus
pHi = 0.5*(lambda0) + (1-lambda0)*(1-normcdf(mu0,mu,sig));

% Choice probability
switch lapseType
    case 'none'
        % No lapse process
        pChoice =  pHi;
    case 'fixed'
        % Fixed proportion of random guesses
        pError = lambda;
        pChoice = pError*0.5 + (1-pError)*pHi;
    case 'inattention'
        % Inattentive guesses according to prior
        pInattention = lambda;
        pChoice = pInattention*(1-normcdf(mu0,muP,sigP)) ...
            + (1-pInattention)*pHi;
end
end


%------------------Simulate choices from ideal observer--------------------
function choices = simulateIdealObserver(params,rules,stim,lapseType)
muP = params.muP;        % Prior mean
sigP = params.sigP;      % Prior standard deviation
sigS = params.sigS;      % Sensory noise
lambda = params.lambda;  % Lapse rate

%Simulate choices
choices = rand(1,length(stim))<pChoiceIdeal([sigS,muP,sigP,lambda], rules, stim, lapseType);
end

%---------Compute mean, S.D. of conditional prior of Xi given Xj=a---------
function [mu, sigma] = getCondPrior(Mu, Sigma, irrDim, a)
switch irrDim
    case 1
        mu = Mu(2) + (Sigma(2,1)/Sigma(1,1))*(a-Mu(1));
        sigma = sqrt(Sigma(2,2)-(Sigma(2,1)^2)/Sigma(1,1));
    case 2
        mu = Mu(1) + (Sigma(1,2)/Sigma(2,2))*(a-Mu(2));
        sigma = sqrt(Sigma(1,1)-(Sigma(1,2)^2)/Sigma(2,2));
end
end

%---------Joint negative log likelihood of data given parameters-----------
function nll = nllJoint(theta, prior, rules, stim, choice, lapseType)
for i = 1:size(stim,2)
    params(1) = theta(2*i-1);
    params(2) = prior{i}.mu;
    params(3) = prior{i}.sig;
    params(4) = theta(2*i);
    
    pChoice = pChoiceIdeal(params,rules{i},stim{i},lapseType);
    thisNLL = -log(choice{i}.*pChoice + ~choice{i}.*(1-pChoice));
    condNLL(i) = sum(thisNLL);
end
nll = sum(condNLL);
end
% 
%-------------------Get empirical prior parameters------------------------

function [jointPrior] = getJointPrior(priorType)
switch priorType
    case 'onset'
        % Onset slope
        jointPrior.Mu = [588.061224489794,38.4702570131941];
        jointPrior.Sigma = [60033.7331794379,-2266.66835833170;-2266.66835833170,589.187367361779];
    case 'offset'
        % Offset slope
        jointPrior.Mu = [588.061224489794,10.0838694529419];
        jointPrior.Sigma = [60033.7331794379,374.221800331715;374.221800331715,88.2037289566602];
    case 'both'
        % % Combined
        jointPrior.Mu = [587.102185380555,0.0100181039239024];
        jointPrior.Sigma = [60108.5861261626,0.0527685193360701;0.0527685193360701,2.14558631529986e-05];
end
end

%-----------------Constrain sigmas to be scaled by sqrt(tHold)-------------
function [c,ceq] = integration_constraint(x)
ceq(1) = x(1)/x(3) - x(7)/x(9);
ceq(2) = x(1)/x(5) - x(7)/x(11);
c = [];
end
