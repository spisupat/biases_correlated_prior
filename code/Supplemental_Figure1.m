
clc
clear
d6005 = load('/data/BurstData/BurstData_600_5.mat');
d60032 = load('/data/BurstData/BurstData_600_32.mat');

d2005 = load('/data/BurstData/BurstData_200_5.mat');
d20032 = load('/data/BurstData/BurstData_200_32.mat');

d200sb5 = load('/data/BurstData/BurstData_200sb_5.mat');
d200sb32 = load('/data/BurstData/BurstData_200sb_32.mat');
load('/data/BurstData/Vcond.mat');

% remove spurious trials that are beyond imposed experimental time constraints
d200sb5.dataset.HoldTime(d200sb5.dataset.HoldTime<=.17 | d200sb5.dataset.HoldTime>=4) = NaN;
d200sb32.dataset.HoldTime(d200sb32.dataset.HoldTime<=.17 | d200sb32.dataset.HoldTime>=4) = NaN;
d2005.dataset.HoldTime(d2005.dataset.HoldTime<=.17 | d2005.dataset.HoldTime>=4) = NaN;
d20032.dataset.HoldTime(d20032.dataset.HoldTime<=.17 | d20032.dataset.HoldTime>=4) = NaN;
d6005.dataset.HoldTime(d6005.dataset.HoldTime<=.57 | d6005.dataset.HoldTime>=4) = NaN;
d60032.dataset.HoldTime(d60032.dataset.HoldTime<=.57 | d60032.dataset.HoldTime>=4) = NaN;



for i=1:length(d200sb5.dataset.Jitter)
    
    h = d200sb5.dataset.HoldTime(i);
    t1 = Vcond(d200sb5.dataset.StimCond(i),d200sb5.dataset.Jitter(i),1);
    t2 = Vcond(d200sb5.dataset.StimCond(i),d200sb5.dataset.Jitter(i),2);
    if h<t2
        bursts(i) = (h-t1)/(t2-t1);
    elseif h>t2
        bursts(i) = 1;
    elseif h==NaN
        bursts(i)=NaN;
    end
end
burstsVec{1} = bursts';


clearvars bursts
for i=1:length(d200sb32.dataset.Jitter)
    
    h = d200sb32.dataset.HoldTime(i);
    t1 = Vcond(d200sb32.dataset.StimCond(i),d200sb32.dataset.Jitter(i),1);
    t2 = Vcond(d200sb32.dataset.StimCond(i),d200sb32.dataset.Jitter(i),2);
    if h<t2
        bursts(i) = (h-t1)/(t2-t1);
    elseif h>t2
        bursts(i) = 1;
    elseif h==NaN
        bursts(i)=NaN;
    end
end
burstsVec{4} = bursts';



clearvars bursts
for i=1:length(d2005.dataset.Jitter)
    
    h = d2005.dataset.HoldTime(i);
    t1 = Vcond(d2005.dataset.StimCond(i),d2005.dataset.Jitter(i),1);
    t2 = Vcond(d2005.dataset.StimCond(i),d2005.dataset.Jitter(i),2);
    t3 = Vcond(d2005.dataset.StimCond(i),d2005.dataset.Jitter(i),3);
    t4 = Vcond(d2005.dataset.StimCond(i),d2005.dataset.Jitter(i),4);
    
    if h<t2
        bursts(i) = (h-t1)/(t2-t1);
    elseif h>t2 && h<t3
        bursts(i) = 1;
    elseif h>t3 &&  h<t4
        bursts(i) =1+ (h-t3)/(t4-t3);
    elseif h>t4 && h<1+t1
        bursts(i) =2;
    elseif h>1+t1
        bursts(i) =2+ (h-(t4+t1))/((t2)-(t1));
    elseif h==NaN
        bursts(i)=NaN;
    end
end
burstsVec{2} = bursts';

clearvars bursts
for i=1:length(d20032.dataset.Jitter)
    
    h = d20032.dataset.HoldTime(i);
    t1 = Vcond(d20032.dataset.StimCond(i),d20032.dataset.Jitter(i),1);
    t2 = Vcond(d20032.dataset.StimCond(i),d20032.dataset.Jitter(i),2);
    t3 = Vcond(d20032.dataset.StimCond(i),d20032.dataset.Jitter(i),3);
    t4 = Vcond(d20032.dataset.StimCond(i),d20032.dataset.Jitter(i),4);
    
    if h<t2
        bursts(i) = (h-t1)/(t2-t1);
    elseif h>t2 && h<t3
        bursts(i) = 1;
    elseif h>t3 &&  h<t4
        bursts(i) =1+ (h-t3)/(t4-t3);
    elseif h>t4 && h<1+t1
        bursts(i) =2;
    elseif h>1+t1
        bursts(i) =2+ (h-(1+t1))/((t2)-(t1));
    elseif h==NaN
        bursts(i)=NaN;
    end
end
burstsVec{5} = bursts';



clearvars bursts
for i=1:length(d6005.dataset.Jitter)
    
    h = d6005.dataset.HoldTime(i);
    t1 = Vcond(d6005.dataset.StimCond(i),d6005.dataset.Jitter(i),1);
    t2 = Vcond(d6005.dataset.StimCond(i),d6005.dataset.Jitter(i),2);
    t3 = Vcond(d6005.dataset.StimCond(i),d6005.dataset.Jitter(i),3);
    t4 = Vcond(d6005.dataset.StimCond(i),d6005.dataset.Jitter(i),4);
    
    if h<t2
        bursts(i) = (h-t1)/(t2-t1);
    elseif h>t2 && h<t3
        bursts(i) = 1;
    elseif h>t3 &&  h<t4
        bursts(i) =1+ (h-t3)/(t4-t3);
    elseif h>t4 && h<1+t1
        bursts(i) =2;
    elseif h>1+t1
        bursts(i) =2+ (h-(1+t1))/((t2)-(t1));
    elseif h==NaN
        bursts(i)=NaN;
    end
end
burstsVec{3} = bursts';

clearvars bursts
for i=1:length(d60032.dataset.Jitter)
    
    h = d60032.dataset.HoldTime(i);
    t1 = Vcond(d60032.dataset.StimCond(i),d60032.dataset.Jitter(i),1);
    t2 = Vcond(d60032.dataset.StimCond(i),d60032.dataset.Jitter(i),2);
    t3 = Vcond(d60032.dataset.StimCond(i),d60032.dataset.Jitter(i),3);
    t4 = Vcond(d60032.dataset.StimCond(i),d60032.dataset.Jitter(i),4);
    
    if h<t2
        bursts(i) = (h-t1)/(t2-t1);
    elseif h>t2 && h<t3
        bursts(i) = 1;
    elseif h>t3 &&  h<t4
        bursts(i) =1+ (h-t3)/(t4-t3);
    elseif h>t4 && h<1+t1
        bursts(i) = 2;
    elseif h>1+t1
        bursts(i) = 2+ (h-(1+t1))/((t2)-(t1));
    elseif h == NaN
        bursts(i)=NaN;
    end
end
burstsVec{6} = bursts';
holdtimeVec = {d200sb5.dataset.HoldTime , d2005.dataset.HoldTime, d6005.dataset.HoldTime, d200sb32.dataset.HoldTime, d20032.dataset.HoldTime, d60032.dataset.HoldTime};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = flip({[182,0,0] ./ 255,[255,7,45] ./ 255,[255,95,146] ./ 255 , [0,0,142] ./ 255, [7,72,255] ./ 255, [95,177,255] ./ 255});

for i = 1:6
    h = holdtimeVec{i};
    b = burstsVec{i};
    
    h = h(b<=2);
    b = b(b<=2);
    
    x = [0:.0001:1.5];
    
    subplot(3,3,i)
    plot(h,b,'.','markeredgecolor',colors{i}); hold on    
    
    break0 = min(h);
    break1 = max(h(b<1));
    break2 = min(h(b>1));
    



    if i ~= 3 && i ~= 6
        mdl1 = fitlm(h(h<break1),b(h<break1),'RobustOpts','cauchy');
        yest1 = mdl1.Coefficients.Estimate(1) + mdl1.Coefficients.Estimate(2).*x;
        y1 = mdl1.Coefficients.Estimate(1) + mdl1.Coefficients.Estimate(2).*x(yest1<1 & x>break0);
        
        plot(x(yest1<1 & x>break0),y1,'k-')
    end
    
    if i ~= 1 && i ~= 4
        mdl2 = fitlm(h(h>break2),b(h>break2),'RobustOpts','cauchy');
        yest2 = mdl2.Coefficients.Estimate(1) + mdl2.Coefficients.Estimate(2).*x;
        y2 = mdl2.Coefficients.Estimate(1) + mdl2.Coefficients.Estimate(2).*x(yest2>1 & yest2<2);
        plot(x(yest2>1 & yest2<2),y2,'k-'); 
        plot(x(yest1>1 & yest2<1),ones(1,length(x(yest1>1 & yest2<1))),'k-')
        plot(x(yest2>2),2.*ones(1,length(x(yest2>2))),'k-')

    end
    

    MedHT(i) = median(h,'omitnan');
    
    if MedHT(i) < x(yest1==max(y1))
        SoundBursts(i) = predict(mdl1,median(h,'omitnan'));
    elseif MedHT(i) > x(yest1==max(y1)) 
        if MedHT(i) < x(yest2==min(y2))
            SoundBursts(i) = 1;
        elseif MedHT(i) > x(yest2==min(y2))
            SoundBursts(i) = predict(mdl2,median(h,'omitnan'));
        end
    end

    ylim([0,2.1])
    xlim([0,1])
    axis square
    box off
ylabel('Hold Time (ms)')
xlabel('Number of Bursts (k)')
end



subplot(3,3,7)

scatter(MedHT,SoundBursts,30,[colors{1};colors{2};colors{3};colors{4};colors{5};colors{6}],'filled')
hold on
plot([0,1],[0,2],'k--')
xlabel('Hold Time (ms)')
xlabel('Number of Bursts (k)')
xlim([0,1])
ylim([0,2])



h = gcf;
saveas(h,'/results/supp_figure_1.pdf')


