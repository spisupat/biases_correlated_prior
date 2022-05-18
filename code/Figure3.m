
clc
clear
load('/data/htdata.mat')
load('/data/htdata2.mat') 
load('/data/paramsValues.mat')
load('/data/dataset_Jane_Pisupati.mat') 


AccMB = [conditionVector{1, 1}(stimulusVector{1,1}~=6 & conditionVector{1, 1}~=2)];



close all
ht = {};
a = {};

colors = {[182,0,0] ./ 255,[255,7,45] ./ 255,[255,95,146] ./ 255 ; [0,0,142] ./ 255, [7,72,255] ./ 255, [95,177,255] ./ 255}
for i=1:3
    for j=1:2
        subplot(4,3,i*j + abs(j-3)*(i-1))
        edges = [.17:.01:1];
        a{j,i} = conditionVector{j,i}(stimulusVector{j,i}~=6 & conditionVector{j,i}~=2);
        histogram(HT{j,i},edges,'facecolor',colors{j,i},'edgecolor','none','normalization','probability')
        hold on
        plot(median(HT{j,i}),.11,'o','color',colors{j,i},'linewidth',1,'markerfacecolor',colors{j,i})
        [IQR]=quantile(HT{j,i},[.25,.75]);
        plot([IQR(1),IQR(2)],[.11,.11],'-','color',colors{j,i},'linewidth',2)
        ylim([0,.12])
        xlim([0,1])
        xlabel('Hold Time (ms)')
    end
end
% close all
% figure
accuracy = [mean(a{1,1}),mean(a{2,1}),mean(a{1,2}),mean(a{2,2}),mean(a{1,3}),mean(a{2,3})];
numTrials = [length(a{1,1}),length(a{2,1}),length(a{1,2}),length(a{2,2}),length(a{1,3}),length(a{2,3})];
medHT = [median(HT{1,1}),median(HT{2,1}),median(HT{1,2}),median(HT{2,2}),median(HT{1,3}),median(HT{2,3})];


subplot(4,3,[10,11])
ciupper = accuracy+1.96.*sqrt(accuracy.*(1-accuracy)./numTrials);
cilower =  accuracy-1.96.*sqrt(accuracy.*(1-accuracy)./numTrials);
plot([medHT;medHT],[cilower;ciupper],'k-','markerfacecolor','k','linewidth',1.5)
hold on
scatter(medHT,accuracy,25,[colors{1,1};colors{2,1};colors{1,2};colors{2,2};colors{1,3};colors{2,3}],'filled')

box off
ylim([.70,.9])
xlim([0,1])



for i = 1:3
    subplot(4,3,i*3)
    xFit = [85:265];
    yFit = PsychoFit(xFit,pV.Fast_Slope(i,:));
    plot(xFit,yFit,'color',colors{1,i})
    hold on
    y = dataset(i+3).data.y;
    x = dataset(i+3).data.x.*100;
    n = dataset(i+3).data.n;
    p = y ./ n;
    sep = sqrt((p.*(1-p))./n);
    U95 = p + 1.96*sep;
    L95 = p - 1.96*sep;
    scatter(x,p,10,colors{1,i},'filled')
    plot([x;x],[U95;L95],'-','color',colors{1,i})

    yFit = PsychoFit(xFit,pV.Slow_Slope(i,:));
    plot(xFit,yFit,'color',colors{2,i})
    y = dataset(i).data.y;
    x = dataset(i).data.x.*100;
    n = dataset(i).data.n;
    p = y ./ n;
    sep = sqrt((p.*(1-p))./n);
    U95 = p + 1.96*sep;
    L95 = p - 1.96*sep;
    scatter(x,p,10,colors{2,i},'filled')
    plot([x;x],[U95;L95],'-','color',colors{2,i})
    xlabel('Duration (ms)','fontsize',8)
    ylabel('Prob. Choose Long','fontsize',8)
    xlim([85,265])
end




h = gcf
saveas(h,'/results/figure_3.pdf')


function y = PsychoFit(x,p)
    mu = p(1);
    sigma = 1./p(2);
    lambda = p(3);
    y = lambda + (1-2*lambda)*normcdf(x,mu,sigma);
end










