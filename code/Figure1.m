function [h,prior] = Figure1()
close all
% clear

% Try loading in summary statistics from current colder
try
    load('/data/Voc_DataMatrix.mat','ex_voc','ex_env','df')
    load('/data/RaterScore.mat'); % load inclusion/exclusion decisions 
    load('/data/SyntheticSounds.mat')
catch
    warning('Vocalization data or Rater Decisions not found in current folder');
end

% clearvars -except datmatrix R df Fs
lowerFreq=14000;  % Lower bound of peak frequency
upperFreq=34000; % Upper bound of peak frequency


% Index frequency bounds
centerFreq = df.PeakFreq;
fbounds = centerFreq>lowerFreq & centerFreq<upperFreq & RaterScores; % Index for vocalizations that are within bounds
dfv = df(fbounds,:); % summary statistics table with non 22kHz vocalizations removed

dur = 1000.*(dfv.thoff-dfv.thon);   % Calculation of duration
onset_slope = (((dfv.Ahon-dfv.Aon)./dfv.AUC)./abs(dfv.thon-dfv.ton)); % Calculation of onset slope
offset_slope = (((dfv.Ahoff-dfv.Aoff)./dfv.AUC)./abs(dfv.toff-dfv.thoff)); % Calculation of offset slope

idx = ~isinf(onset_slope) & ~isinf(offset_slope);
onset_slope = onset_slope(idx);
offset_slope = offset_slope(idx);
dur = dur(idx);

AUC_slow_condition = trapz([1:Synthetic.Fs]./Synthetic.Fs,Synthetic.Slow)/2;
AUC_fast_condition = trapz([1:Synthetic.Fs]./Synthetic.Fs,Synthetic.Fast)/2;
slowDerivative =(diff(Synthetic.Slow)./AUC_slow_condition).*Synthetic.Fs; % Derivative of synthetic sounds (Fast slope condition)
fastDerivative = (diff(Synthetic.Fast)./AUC_fast_condition).*Synthetic.Fs; % Derivative of synthetic sounds (Slow slope condition)
 
x = [1:Synthetic.Fs]./Synthetic.Fs;

% Onset/offset slopes of synthetic sounds (both slope conditions: fast and slow)
slow = max(slowDerivative); % Slope rate at inflection point of synthetic sounds
fast = max(fastDerivative);    % Slope rate at inflection point of synthetic sounds


given = [slow, fast]; % Conditional slope values from synthetic sounds used from the discrimination task
color = {'r','b'};  % Color scheme (slope conditions: red=slow; blue=fast)
slopes = [onset_slope, offset_slope]; % Matrix for all three slope variables

model_type = {'onset_slope', 'offset_slope', 'joint_slope'}'; % Slope variables used for different prior types
parameter_label = {'mu','sigma2','mu|slow','sigma2|slow','mu|fast','sigma2|fast'}; % Labels of parameters for bivariate and conditional gaussians for empirical prior (mu = means, sigma2 = variance)
prior = array2table(NaN(length(model_type),length(parameter_label)),'variablenames',parameter_label,'rownames',model_type); % construct table

prior.('mu') = {mean([dur,onset_slope],'omitnan'),mean([dur,offset_slope],'omitnan'),mean([dur,onset_slope,offset_slope],'omitnan')}'; % Means of multivariate gaussians
prior.('sigma2') = {cov([dur,onset_slope],'omitrows'),cov([dur,offset_slope],'omitrows'), cov([dur,onset_slope,offset_slope],'omitrows')}'; % Covariances of multivariate gaussians
for mdl_idx =1:length(model_type)
    [prior.('mu|slow')(mdl_idx), prior.('sigma2|slow')(mdl_idx)] = condGauss(prior.('mu'){mdl_idx},prior.('sigma2'){mdl_idx},slow);  % Conditional Gaussian parameters for slow slope condition 
    [prior.('mu|fast')(mdl_idx), prior.('sigma2|fast')(mdl_idx)] = condGauss(prior.('mu'){mdl_idx},prior.('sigma2'){mdl_idx},fast); %Conditional Gaussian parameters for fast slope condition 
end
 

%%%%%%%%%%%%%%%%%%%%%%%
%  Large figure panel %
%%%%%%%%%%%%%%%%%%%%%%%

figure('units','normalized','outerposition',[0 0 1 1])
CondLab = {'Onset','Offset'};
x = [0:1500]';

% Overlay regression line w/ confidence intervals on scatter plot

    x_domain =  [0:5:1500];
    y_domain =  [0:2.5:225];
for i = 1:2
    subplot(4,4,8+i)
    [yhat,ci] = predict(fitlm(dur,slopes(:,i)),x,'Alpha',0.05);
    h = scatter(dur,slopes(:,i),40,'.','filled');
    hold on
    h.MarkerFaceColor = [.6,.6,.6];
    h.MarkerEdgeColor = [.6,.6,.6];
    plot(x,yhat,'k-','linewidth',2)
    plot(x,ci,'k:','linewidth',1)
    ylim([0,225])
    xlim([0,1500])
    set(gca,'TickDir','out');
    colorbar
    set(colorbar,'visible','off')  % colorbar adjustment to align subplots
    xlabel('Plateau Duration (ms)')
    ylabel('Slope (A/s)')
    title(CondLab{i});
    axis square

    % Plot bivariate probability density functions as a smoothed surface
    ax=subplot(4,4,12+i);
    mu = prior.('mu'){i}; sigma = prior.('sigma2'){i};
    [X, Y] = meshgrid(x_domain,y_domain);
    variable_space = [X(:), Y(:)];
    Px = mvnpdf(variable_space,mu,sigma);
    Px = reshape(Px,length(y_domain),length(x_domain));
    Px = Px./max(max(Px));
    s = surf(x_domain,y_domain,zeros(size(Px)),Px);
    s.EdgeColor = 'none';
    hold on
    set(gca,'YDir','normal')
    view(2)
    axis square
    colormap(ax,[[1:(-1/256):0]',[1:(-1/256):0]',[1:(-1/256):0]'])
    set(gca,'TickDir','out');
    ylim([0,225])
    ylabel('Slope (A/ms)')
    xlabel('Duration (ms)')
    title(CondLab{i});
end

% Plot conditional gaussians for each type of prior
for i = [1:3]
subplot(4,4,[11,12,15,16])
x_domain =  [-2500:2:2500];
y = normpdf(x_domain,prior.('mu|slow')(i)-175,sqrt(prior.('sigma2|slow')(i))); hold on
plot(x_domain,(y./max(y))./1.3+i,'r-','linewidth',2)
plot([prior.('mu|slow')(i)-175, prior.('mu|slow')(i)-175],[i, max((y./max(y))./1.3+i)],'r--','linewidth',1.5)

y = normpdf(x_domain,prior.('mu|fast')(i)-175,sqrt(prior.('sigma2|fast')(i)));
plot(x_domain,(y./max(y))./1.3+i,'b-','linewidth',2)
plot([prior.('mu|fast')(i)-175, prior.('mu|fast')(i)-175],[i, max((y./max(y))./1.3+i)],'b--','linewidth',1.5)

plot([0,3000],[i,i],'-','color',[.6,.6,.6])
title('Priors')
end
yticks([1:3])
ylim([0.75,4])
xlim([-2500,2500])
yticklabels({'Onset','Offset','Both(Joint)'})
xlabel('Duration Relative to Categorical Boundary')

%Plot spectrogram of vocalization
ax = subplot(4,5,[6,7]);
xf = 1:length(ex_voc);
Fs = 250000;
Time = xf/Fs;
spectrogram(ex_voc,128*15,120*15,128*15,Fs,'yaxis') % Spectrogram plot display
hold on
set(colorbar,'visible','on')  % Display color bar
colormap(ax,[[1:(-1/256):0]',[1:(-1/256):0]',[1:(-1/256):0]'])
xticks([0:3])
xlim([0 3])
yticks([])
ylim([0 80])
xlabel('Time (s)','Fontsize',12)
ylabel('Frequency (kHz)','fontsize',12)
caxis([-110 -90])        % set colorbar axis for contrast
box off
yticks([0:10:80])

%Plot raw signal of example vocalizations
subplot(4,5,[1,2])
plot(Time,ex_voc,'k-')
colorbar
set(colorbar,'visible','off')  % colorbar aligns spectrogram with signal
ylabel('Amplitude (a.u.)','fontsize',12)
xlabel('Time (s)','Fontsize',12)
xticks([0:1:3])
yticks([-1:1])
box off

% plot envelope of first example vocalization
a = 834.5;
Idx = 500;
subplot(4,4,[3,4])
max_env = max(ex_env(Time<1.75 & Time>.2));
plot(Time(Time<1.75 & Time>.2)-.2,ex_env(Time<1.75 & Time>.2)'./max_env,'-','linewidth',1.5,'color',[.2 .6 .6])
hold on
xlabel('Time (s)','Fontsize',12)

%plot onset, offset and half-max points onto envelope
tonex = df.ton(Idx)-a;
Aonex = df.Aon(Idx)./max_env;
thonex = df.thon(Idx)-a;
Ahonex = df.Ahon(Idx)./max_env;

toffex = df.toff(Idx)-a;
Aoffex = df.Aoff(Idx)./max_env;
thoffex = df.thoff(Idx)-a;
Ahoffex = df.Ahoff(Idx)./max_env;
for call=1:length(Idx)
plot([tonex(call),thonex(call),thoffex(call),toffex(call)],[Aonex(call),Ahonex(call),Ahoffex(call),Aoffex(call)],'ko-','linewidth',2)
end
xlim([0,1.5])
box off

% plot envelope of synthetic sounds
subplot(4,4,[7,8])
Xs = [1:length(Synthetic.Short)]./Synthetic.Fs;
plot(Xs(Xs<.5),Synthetic.Short(Xs<.5),'-','linewidth',2,'color',[.3 .7 .7]); hold on
plot(Xs(Xs<.5),Synthetic.Middle(Xs<.5),'-','linewidth',2,'color',[.2 .6 .6])
plot(Xs(Xs<.5),Synthetic.Long(Xs<.5),'-','linewidth',2,'color',[.1 .5 .5])
xlim([0,1.5])
ylim([0,1.1])
xlabel('Time (s)','Fontsize',12)
box off

clearvars -except prior
h=gcf;
set(h,'PaperSize',[20 10]); %set the paper size to what you want  
%print(h,'/results/figure_panel.pdf') % then print it

saveas(h,'/results/figure_1.pdf')
end 


function [condMu,condSigma] = condGauss(Mu,Sigma,a)

% conditional gaussian calculation (handles 2-3 variables)
if length(Mu) == 2
    condMu = Mu(1)+Sigma(1,2)/Sigma(2,2)*(a-Mu(2));
    condSigma = Sigma(1,1)-(Sigma(1,2)/Sigma(2,2))*Sigma(2,1);
elseif length(Mu) == 3
    
    Mu1 = Mu(1:2)+(Sigma(3,1:2)/Sigma(3,3))*(a-Mu(3));
    Sigma1 = Sigma(1:2,1:2)-(Sigma(1:2,3)/Sigma(3,3))*Sigma(3,1:2);
    
    condMu = Mu1(1)+(Sigma1(1,2)/Sigma1(2,2))*(a-Mu1(2));
    condSigma = Sigma1(1,1)-(Sigma1(1,2)/Sigma1(2,2))*Sigma1(2,1);
    
end

end
