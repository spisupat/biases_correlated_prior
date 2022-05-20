
load('/data/BootstrapDataMat.mat')

addpath('/code/Functions/Violin.m')
addpath('/code/Functions/violinplot.m')
x=[50:.5:300];


for i = 1:3
     if i == 1
        k=2;
    elseif i == 2
        k=3;
     elseif i == 3
        k=1;
     end
    
    Yl5(:,i) = [paramsSim(:,i,3)+paramsSim(:,i,4)];
    
    Yl32(:,i) = [paramsSim(:,i+3,3)+paramsSim(:,i+3,4)];
    
    Ym5(:,i) = [paramsSim(:,i,1)];
    
    Ym32(:,i) = [paramsSim(:,i+3,1)];
    
    Ys5(:,i) = [1./paramsSim(:,i,2)];
    
    Ys32(:,i) = [1./paramsSim(:,i+3,2)];
    X1(:,i) = k.*ones(400,1);
end




Y = Ym5(:);
X = X1(:);
subplot(3,2,[1])
violinplot(Y,X,'ViolinColor',[.7,.7,.7],'violinalpha',.01,'BoxColor',[1,0,0],'EdgeColor',[1,0,0]);
xticks([1:3])
xticklabels({'200ms SB','200ms MB','600ms MB'})
xtickangle(45)
ylim([160,205])
yticks([160:15:205])
ylabel('Bias (\mu)','fontsize',18)
xlabel('Average Bursts','fontsize',18)
box off

Y = Ym32(:);
subplot(3,2,[2])
violinplot(Y,X,'ViolinColor',[.7,.7,.7],'violinalpha',.01,'BoxColor',[0,0,1],'EdgeColor',[0,0,1]);
xticks([1:3])
xticklabels({'200ms SB','200ms MB','600ms MB'})
xtickangle(45)
ylim([160,205])
yticks([160:15:205])
box off
%
Y = Ys5(:);
X = X1(:);
subplot(3,2,[3])
violinplot(Y,X,'ViolinColor',[.7,.7,.7],'violinalpha',.01,'BoxColor',[1,0,0],'EdgeColor',[1,0,0]);
xticks([1:3])
xticklabels({'200ms SB','200ms MB','600ms MB'})
xtickangle(45)
ylim([30,90])
yticks([30:10:90])
ylabel('Inverse Sensitivity (\sigma)','fontsize',18)
box off

Y = Ys32(:);
subplot(3,2,[4])
violinplot(Y,X,'ViolinColor',[.7,.7,.7],'violinalpha',.01,'BoxColor',[0,0,1],'EdgeColor',[0,0,1]);
xticks([1:3])
xticklabels({'200ms SB','200ms MB','600ms MB'})
xtickangle(45)
ylim([30,90])
yticks([30:10:90])
box off
%

Y = Yl5(:);
X = X1(:);
subplot(3,2,[5])
violinplot(Y,X,'ViolinColor',[.7,.7,.7],'violinalpha',.01,'BoxColor',[1,0,0],'EdgeColor',[1,0,0]);
xticks([1:3])
xticklabels({'200ms SB','200ms MB','600ms MB'})
xtickangle(45)
ylim([0,.2])
yticks([0:.05:.2])
ylabel('Lapse Rate (\lambda+\gamma)','fontsize',18)
box off

Y = Yl32(:);
subplot(3,2,[6])
violinplot(Y,X,'ViolinColor',[.7,.7,.7],'violinalpha',.01,'BoxColor',[0,0,1],'EdgeColor',[0,0,1]);
xticks([1:3])
xticklabels({'200ms SB','200ms MB','600ms MB'})
xtickangle(45)
ylim([0,.2])
yticks([0:.05:.2])
box off

h = gcf;
saveas(h,'/results/figure_4.pdf')

