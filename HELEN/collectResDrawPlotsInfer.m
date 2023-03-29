%% 
% load basic data
clear;

continentsAll = {'Europe','Europe','North_America','North_America','Asia','Asia','Africa','Africa','Oceania','South_America','South_America','Europe','Europe','Europe','Europe','Europe'};
countriesAll = {'United_Kingdom','Germany','USA','Canada','Japan','India','South Africa','Kenya','Australia','Brazil','Peru','France','Sweden','Spain','Denmark','Italy'};
countriesAllPrint = countriesAll;
countriesAllPrint{1} = 'UK';
countriesAllPrint{7} = 'RSA';

addpath 'C:\gurobi811\win64\matlab'

mut_alpha = [501 570 614 681 716 982 1118];
mut_beta = [80 215 417 484 501 614 701];
mut_gamma = [18 20 26 138 190 417 484 501 614 655 1027 1176];
mut_delta = [19 142 158 452 478 614 681 950];
mut_omicron = [67 95 145 212 339 371 373 375 417 440 446 477 478 484 493 496 498 501 505 547 614 655 679 681 764 796 856 954 969 981];
mut_lambda = [75 76 452 490 614 859];
mut_mu = [95 145 346 484 501 614 681 950];
mut_theta = [265 484 501 614 681 1092 1101 1176];
mut_eta = [52 67 484 614 677 888];
mut_kappa = [142 154 452 484 614 681 1071];
mut_zeta = [484 614 1176];
mut_epsilon = [13 152 452 614];
mut_iota = [5 95 614];

% nMut = size(M,2);
nMut = 1273;
voc = false(13,nMut);
voc(1,mut_alpha) = 1;
voc(2,mut_beta) = 1;
voc(3,mut_gamma) = 1;
voc(4,mut_delta) = 1;
voc(5,mut_omicron) = 1;
voc(6,mut_lambda) = 1;
voc(7,mut_mu) = 1;
voc(8,mut_theta) = 1;
voc(9,mut_eta) = 1;
voc(10,mut_kappa) = 1;
voc(11,mut_zeta) = 1;
voc(12,mut_epsilon) = 1;
voc(13,mut_iota) = 1;

voc_design = cell(1,13);
voc_design{1} = datetime('12-18-2020','InputFormat','MM-dd-yyyy');
voc_design{2} = datetime('12-18-2020','InputFormat','MM-dd-yyyy');
voc_design{3} = datetime('01-11-2021','InputFormat','MM-dd-yyyy');
voc_design{4} = datetime('05-11-2021','InputFormat','MM-dd-yyyy');
voc_design{5} = datetime('11-26-2021','InputFormat','MM-dd-yyyy');
voc_design{6} = datetime('07-14-2021','InputFormat','MM-dd-yyyy');
voc_design{7} = datetime('08-30-2021','InputFormat','MM-dd-yyyy');
voc_design{8} = datetime('03-24-2021','InputFormat','MM-dd-yyyy');
voc_design{9} = datetime('03-17-2021','InputFormat','MM-dd-yyyy');
voc_design{10} = datetime('04-04-2021','InputFormat','MM-dd-yyyy');
voc_design{11} = datetime('03-17-2021','InputFormat','MM-dd-yyyy');
voc_design{12} = datetime('03-05-2021','InputFormat','MM-dd-yyyy');
voc_design{13} = datetime('03-24-2021','InputFormat','MM-dd-yyyy');

voc_design_fig = cell(1,13);
voc_design_fig{1} = datetime('12-18-2020','InputFormat','MM-dd-yyyy');
voc_design_fig{2} = datetime('12-18-2020','InputFormat','MM-dd-yyyy');
voc_design_fig{3} = datetime('01-11-2021','InputFormat','MM-dd-yyyy');
voc_design_fig{4} = datetime('05-11-2021','InputFormat','MM-dd-yyyy')+3;
voc_design_fig{5} = datetime('11-26-2021','InputFormat','MM-dd-yyyy');
voc_design_fig{6} = datetime('07-14-2021','InputFormat','MM-dd-yyyy');
voc_design_fig{7} = datetime('08-30-2021','InputFormat','MM-dd-yyyy');
voc_design_fig{8} = datetime('03-24-2021','InputFormat','MM-dd-yyyy');
voc_design_fig{9} = datetime('03-17-2021','InputFormat','MM-dd-yyyy');
voc_design_fig{10} = datetime('04-04-2021','InputFormat','MM-dd-yyyy');
voc_design_fig{11} = datetime('03-17-2021','InputFormat','MM-dd-yyyy');
voc_design_fig{12} = datetime('03-05-2021','InputFormat','MM-dd-yyyy');
voc_design_fig{13} = datetime('03-24-2021','InputFormat','MM-dd-yyyy');

voc_colors = [1 0 0; 0 0 1; 0.5 0.5 0;  0.5 0 0.5; 0.1 1 0.2; 0 0 0; 0.3010 0.7450 0.9330; 0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];
voc_names = {'alpha','beta','gamma', 'delta','omicron','lambda','mu','theta','eta','kappa','zeta','epsilon','iota'};


time_test_begin = datetime('05-01-2020','InputFormat','MM-dd-yyyy');
time_test_end = datetime('11-01-2021','InputFormat','MM-dd-yyyy');
time_test_end_extend = datetime('04-22-2022','InputFormat','MM-dd-yyyy');
time_points = time_test_begin:caldays(15):time_test_end;
time_points_extend = time_test_begin:caldays(15):time_test_end_extend;

load(['HELEN data' filesep 'Preprocessed input data' filesep 'data_all_2022.mat']);
mutCounts = sum(M_all,2);
%% 
%analysis of inferred haplotypes; f-score plots for countries separately
eps = 0; %0.001;
fThr = 0 % 0.5;
interv = 1 % 3;
mult = 10;
scatterPointsT = [];
scatterPointsF = [];
scatterPointsV = [];
scatterPointsFreq = [];
predictTime = -Inf*ones(length(countriesAll),size(voc,1)-3);
predictFscore = -Inf*ones(length(countriesAll),size(voc,1)-3);
predictTime1 = -Inf*ones(length(countriesAll),size(voc,1)-3);
predictFreq = -Inf*ones(length(countriesAll),size(voc,1)-3);
predictFreqVocInfer = -Inf*ones(length(countriesAll),size(voc,1)-3);
predictPrevVocInfer = -Inf*ones(length(countriesAll),size(voc,1)-3);
matchDistAll = zeros(length(countriesAll),length(time_points));
matchDistAllExtend = zeros(length(countriesAll),length(time_points));
numEdgesAll = zeros(length(countriesAll),length(time_points));
growDetected = false(1,size(M_all,1));

for c = 1:length(countriesAll)
%     c = 16
    country = countriesAll{c}
    load(['HELEN data' filesep 'Inferred variants' filesep countriesAll{c} '_densestFixSizeExt.mat']);
    load(['HELEN data' filesep 'Preprocessed input data' filesep 'vocs_' countriesAll{c} '.mat']);
    
    fData = NaN*ones(length(time_points),size(voc,1)-3);
    sensData = NaN*ones(length(time_points),size(voc,1)-3);
    matchDist = zeros(1,length(time_points));
    netDens = zeros(1,length(time_points));
    predictVocDynam = cell(1,length(time_points));
    for t = 1:length(time_points)
        t

        [valMaxF,indMaxF] = max(stats{t}.allF,[],2);
        [valMaxSens,indMaxSens] = max(stats{t}.allSens,[],2);
        for v = unique(indMaxF)'
            ind = find((indMaxF'==v)&(stats{t}.vocFreq >= eps));
            if ~isempty(ind)
                fData(t,v) = max(valMaxF(ind));
            end
        end
        for v = unique(indMaxSens)'
            ind = find((indMaxSens'==v)&(stats{t}.vocFreq >= eps));
            if ~isempty(ind)
                sensData(t,v) = max(valMaxSens(ind));
            end
        end
        if ~isempty(stats{t}.vocInfer)
            matchDist(t) = stats{t}.matchDist;
        else
            matchDist(t) = -1;
        end
        if ~isempty(stats{t}.allF)
            for s = 1:10
                if (max(stats{t}.allF(:,s))>=0.8) && (isinf(predictTime(c,s)))
                    predictTime(c,s) = -days(time_points(t)-voc_design{s});
                    predictFreq(c,s) = freq_voc_all(s,t);
                    predictFscore(c,s) = fData(t,s);
                    mv = find(stats{t}.allF(:,s)==max(stats{t}.allF(:,s)),1,'first');
                    predictFreqVocInfer(c,s) = stats{t}.vocInferCummFreq(mv);
                    predictPrevVocInfer(c,s) = stats{t}.vocInferPrev(mv);
                end
                if (max(stats{t}.allF(:,s))>=0.8) && (isinf(predictTime1(c,s)))
                    tm1perc = time_points(min([find(relPrev_voc_all(s,:)>=0.01,1,'first') length(time_points)]));
                    if tm1perc == time_points(end)
                        tm1perc = max([tm1perc voc_design{s}]);
                    end
                    predictTime1(c,s) = -days(time_points(t)-tm1perc);
                end
            end
        end
        
%         if ~isempty(stats{t}.allF)
%             for v = 1:length(stats{t}.vocInfer)
%                 if valMaxF(v)>=0.8
%                     matchDistAllExtend(c,t) = matchDistAllExtend(c,t) + stats{t}.vocFreq(v)*valMaxF(v);
%                 else
%                     mutCountsV = sum(M_all(:,stats{t}.vocInfer{v}),2);
%                     rec = mutCountsV./mutCounts;
%                     prec = mutCountsV/length(stats{t}.vocInfer{v});
%                     fscore = 2*(prec.*rec)./(prec+rec);
%                     ind = (prec >= 0.8)';
%                     ind_future = ind&(sampTimeSeqD_all > time_points(t));
%                     ind_past = ind&(sampTimeSeqD_all <= time_points(t));
%                     ind_detected = ind&growDetected;
%                     if ((sum(ind_future)>=mult*sum(ind_past))&&(sum(ind_future)>=mult))||(sum(ind_detected)>0)
%                         matchDistAllExtend(c,t) = matchDistAllExtend(c,t) + stats{t}.vocFreq(v)*max(fscore(ind_future|ind_detected));
%                         growDetected(ind_future) = 1;
%                     else
%                         matchDistAllExtend(c,t) = matchDistAllExtend(c,t) + stats{t}.vocFreq(v)*valMaxF(v);
%                     end
%                 end
%             end
%         else
%             matchDistAllExtend(c,t) = -1;
%         end
        scatterPointsT = [scatterPointsT repmat(time_points(t),1,sum(~isnan(fData(t,:))))];
        scatterPointsF = [scatterPointsF fData(t,~isnan(fData(t,:)))];
        scatterPointsV = [scatterPointsV find(~isnan(fData(t,:)))];
    end
    
    
    Data = fData;
    Data(isnan(Data)) = 0;
    Data(Data<=fThr) = 0;
    matchDistAll(c,:) = matchDist;
    

    R = ceil(length(time_points)/interv);
    indT = mat2cell((1:interv*R)',interv*ones(1,R));
    indT{end} = indT{end}(indT{end}<=length(time_points));
    time_points_adj = NaT(1,length(indT));
    Data_adj = zeros(length(indT),size(Data,2));
    for i = 1:length(indT)
        time_points_adj(i) = mean(time_points(indT{i}));
        Data_adj(i,:) = max(Data(indT{i},:),[],1);
    end
    [Data_sorted,Data_var] = sort(Data_adj,2,'descend'); 
    Data_sorted = Data_sorted(:,sum(Data_sorted,1)>0);
    Data_var = Data_var(:,sum(Data_sorted,1)>0);
    
    figure
    b = bar(time_points_adj,Data_sorted,'BarWidth', 1,'FaceColor','flat','EdgeColor','k');
    for i = 1:length(b)
        for j = 1:size(Data_sorted,1)
            b(i).CData(j,:) = voc_colors(Data_var(j,i),:);
        end
    end
    title(strrep(countriesAllPrint{c},'_',' '))
    datetick('x',28,'keepticks');
    xlim([time_points_adj(1)-days(15) voc_design{5}+days(15)])
    grid on
    ylabel('f-score','fontweight','bold','fontsize',14)
    ylim([0 1])
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'FontWeight','bold','FontSize',14)
    set(gcf, 'Position',  [100, 100, 1400, 500])


    vocs_found = find(sum(Data_adj,1)>0);
    for s = vocs_found
        xl = xline(voc_design_fig{s},'--','Color',voc_colors(s,:),'LineWidth',3);
    %     xl.LabelOrientation = 'horizontal';
    end
    hold on
    voc_leg = [];
    for i = 1:length(vocs_found)
        voc_leg(i) = plot(NaN,NaN,'Color',voc_colors(vocs_found(i),:),'LineWidth',4);
    end
    hold off
    legend(voc_leg,voc_names(vocs_found),'location', 'Northeastoutside')

%     exportgraphics(gcf,['results' filesep 'figures' filesep countriesAll{c} '_haplotypes.png'],'Resolution',600)

    t1 = find(matchDist~=-1,1,'first');
    median(matchDist(t1:end))    
end

val = 0:0.05:1;
ppv = zeros(1,length(val));
for i = 1:length(val)
    ppv(i) = sum(matchDistAll(:)>=val(i))/sum(matchDistAll(:)~=-1);
end
mean(matchDistAll(matchDistAll~=-1))

val = 0:0.05:1;
ppv = zeros(1,length(val));
for i = 1:length(val)
    ppv(i) = sum(matchDistAllExtend(:)>=val(i))/sum(matchDistAllExtend(:)~=-1);
end

mean(matchDistAllExtend(matchDistAllExtend~=-1))

predictTime_voc = predictTime(:,1:5);
predictTime1_voc = predictTime1(:,1:5);
predictTime_voi = predictTime(:,6:10);
predictTime1_voi = predictTime1(:,6:10);
median(predictTime_voc(predictTime_voc>0))
median(predictTime1_voc(predictTime1_voc > 0))
median(predictTime_voi(predictTime_voi>0))
median(predictTime1_voi(predictTime1_voi > 0))
matchDistAllExtend1 = matchDistAllExtend;
matchDistAllExtend1(matchDistAllExtend1 == - 1) = nan;
%% 
% plot specificity
y = matchDistAllExtend(:);
x = repmat(cellstr(time_points),length(countriesAll),1);
x = x(:);
x(y==-1) = [];
y(y==-1) = [];
figure
set(gcf, 'Position',  [100, 100, 1800, 450])
boxplot(y,x,'PlotStyle','compact','LabelOrientation','horizontal')
ylabel('matching similarity','fontweight','bold','fontsize',12)
ylim([0 1])
grid on
hold on
med = zeros(1,length(time_points));
for i = 1:length(time_points)
    med(i) = median(y(x==time_points(i)));
end
plot(med,'r-','LineWidth',2)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'FontWeight','bold')
labels = cellstr(time_points);
labels(2:2:end) = {' '}; % remove every other one
ax = gca;
ax.XAxis.TickLabels = labels; % set
xl = xlim;
yl = ylim;
text(3.5*xl(1),0.99*yl(2),'(f)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
% exportgraphics(gca,['results' filesep 'figures' filesep 'specificity_infer.png'],'Resolution',600) 

%% 
% f-scores of inferred haplotypes bar plot
ft = 11;
nVoc = 10;
fig = figure
tile = tiledlayout(nVoc,1);
axes = [];
for v = 1:nVoc
    data = zeros(length(time_points),length(countriesAll));
    for t = 1:length(time_points)
        ind = find(scatterPointsT==time_points(t)&(scatterPointsV==v));
        [data(t,1:length(ind)),indSort] = sort(scatterPointsF(ind),'descend');
    end
    k = find(sum(data,1)>0,1,'last');
    data = data(:,1:k);
%     figure
    axes(v) = nexttile;
    b = bar(axes(v),time_points,data,'BarWidth', 1,'FaceColor','flat','EdgeColor','k','LineWidth',0.25);
    grid on
    for i = 1:length(b)
        for j = 1:size(data,1)
            b(i).CData(j,:) = voc_colors(v,:);
        end
    end
%     ylabel('f-score','fontweight','bold','fontsize',14)
%     ylim([0 1])
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'FontWeight','bold','FontSize',ft)
    set(gcf, 'Position',  [50, 50, 1900, 1200])
    xl = xline(voc_design_fig{v},'--','Color','k','LineWidth',3);
    yl = yline(0.8,':','Color','k','LineWidth',0.5);
    if v ~= nVoc
        xticklabels(axes(v),{});
    end
    if v == nVoc
        hold on
        voc_leg = [];
        for i = 1:nVoc
            voc_leg(i) = plot(NaN,NaN,'Color',voc_colors(i,:),'LineWidth',4);
        end
        hold off
        legend(voc_leg,voc_names,'location', 'Southwest','FontSize',ft)
    end
    xl = xlim;
    yl = ylim;
    if v == 1
        text(voc_design{5},0.99*yl(2),'(a)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top');
    end
    yticks([0 0.5 0.8]);
end
linkaxes(axes,'x');
% title(tile,'My Title')
% xlabel(tile,'time')
ylabel(tile,'fscore','fontweight','bold','fontsize',ft)
xl = xlim;
xlim([xl(1) voc_design{5} + days(7)]) 
tile.TileSpacing = 'none';
% exportgraphics(fig,['results' filesep 'figures' filesep 'inferFscoreAll.png'],'Resolution',600)


%% 
% inferred communities forecasting depths and frequency bar plot
ft = 11;
vocFound = find(sum(~isinf(predictTime),1)>0);
w0 = 4;

nRow = 1;
nCol = 4;
figure
[ha, pos] = tight_subplot(nRow,nCol,[.045 .03],[.1 0.05],[.1 .1]);
set(gcf, 'Position',  [100, 100, 2500, 570])

% figure
axes(ha(1));
title('(a) Forecasting depth wrt 1% prevalence time')
% set(gcf, 'Position',  [100, 100, 1190, 805])
xNext = 1;
tks = zeros(1,length(vocFound));
for v = 1:length(vocFound)
%     ptSorted = sort(predictTime(~isinf(predictTime(:,v)),v),'descend');
    ptSorted = sort(predictTime1(:,vocFound(v)),'descend');
    nPred = sum(~isinf(ptSorted));
    if ~isempty(ptSorted)
        b = bar(xNext,ptSorted,1);
        grid on
        for i = 1:length(b)
%             b(i).CData = voc_colors(v,:);
            b(i).FaceColor = voc_colors(vocFound(v),:);
            if ptSorted(i) == 0
                hold on
                bw = b(i+1).XEndPoints - b(i).XEndPoints;
                rectangle('Position',[b(i).XEndPoints - bw/2,-w0,bw,w0],'FaceColor',voc_colors(vocFound(v),:));
            end
        end
    end
    xNext = b(nPred).XEndPoints + 0.5;
    tks(v) = (b(1).XEndPoints + b(nPred).XEndPoints)/2;
    hold on
end
xl = xlim;
xlim([xl(1) b(nPred).XEndPoints+0.1])
xticks(tks)
xticklabels(voc_names(vocFound))
xtickangle(45)
ax = gca;
ax.FontSize = ft;
ax.FontWeight = 'bold';
ylabel('FD_c^{prev}, days')
yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = 0.25 + p(1); 
set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(b)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')

% figure
axes(ha(2));
title('(b) Forecasting depth wrt designation time')
% set(gcf, 'Position',  [100, 100, 1190, 805])
xNext = 1;
tks = zeros(1,length(vocFound));
for v = 1:length(vocFound)
    ptSorted = sort(predictTime(:,vocFound(v)),'descend');
    nPred = sum(~isinf(ptSorted));
    if ~isempty(ptSorted)
        b = bar(xNext,ptSorted,1);
        grid on
        for i = 1:length(b)
%             b(i).CData = voc_colors(v,:);
            b(i).FaceColor = voc_colors(vocFound(v),:);
            if ptSorted(i) == 0
                hold on
                bw = b(i+1).XEndPoints - b(i).XEndPoints;
                rectangle('Position',[b(i).XEndPoints - bw/2,-w0,bw,w0],'FaceColor',voc_colors(vocFound(v),:));
            end
        end
    end
    xNext = b(nPred).XEndPoints + 0.5;
    tks(v) = (b(1).XEndPoints + b(nPred).XEndPoints)/2;
    hold on
end
xl = xlim;
xlim([xl(1) b(nPred).XEndPoints+0.1])
xticks(tks)
xticklabels(voc_names(vocFound))
xtickangle(45)
ax = gca;
ax.FontSize = ft;
ax.FontWeight = 'bold';
ylabel('FD_c^{des}, days')
ylabel('FD_c^{des}, days')
yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = 0.2 + p(1); 
set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(c)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')


% figure
axes(ha(3));
title('(c) Cummulative frequencies at times of variant calling')
% set(gcf, 'Position',  [100, 100, 1190, 805])
xNext = 1;
tks = zeros(1,length(vocFound));
for v = 1:length(vocFound)
    ptSorted = sort(predictFreqVocInfer(:,vocFound(v)),'descend');
    nPred = sum(~isinf(ptSorted));
    if ~isempty(ptSorted)
        b = bar(xNext,ptSorted,1);
        grid on
        for i = 1:length(b)
%             b(i).CData = voc_colors(v,:);
            b(i).FaceColor = voc_colors(vocFound(v),:);
            if ptSorted(i) == 0
                hold on
                p = plot([b(i).XEndPoints - bw/2  b(i).XEndPoints+bw/2],[0.000001 0.000001],':');
                p.LineWidth = 5;
                p.Color = voc_colors(vocFound(v),:);
            end
        end
    end
    xNext = b(nPred).XEndPoints + 0.5;
    tks(v) = (b(1).XEndPoints + b(nPred).XEndPoints)/2;
    hold on
end
xl = xlim;
xlim([xl(1) b(nPred).XEndPoints+0.1])
xticks(tks)
xticklabels(voc_names(vocFound))
xtickangle(45)
ax = gca;
set(gca,'yscale','log')
ax.FontSize = ft;
ax.FontWeight = 'bold';
ylabel('cummulative frequency (log scale)')
yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = 0.05 + p(1); 
set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(d)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')

% figure
axes(ha(4));
title('(d) Prevalences at times of variant calling')
% set(gcf, 'Position',  [100, 100, 1190, 805])
xNext = 1;
tks = zeros(1,length(vocFound));
for v = 1:length(vocFound)
    ptSorted = sort(predictPrevVocInfer(:,vocFound(v)),'descend');
    nPred = sum(~isinf(ptSorted));
    if ~isempty(ptSorted)
        b = bar(xNext,ptSorted,1);
        grid on
        for i = 1:length(b)
%             b(i).CData = voc_colors(v,:);
            b(i).FaceColor = voc_colors(vocFound(v),:);
            if ptSorted(i) == 0
                hold on
                p = plot([b(i).XEndPoints - bw/2  b(i).XEndPoints+bw/2],[0.000001 0.000001],':');
                p.LineWidth = 5;
                p.Color = voc_colors(vocFound(v),:);
            end
        end
    end
    xNext = b(nPred).XEndPoints + 0.5;
    tks(v) = (b(1).XEndPoints + b(nPred).XEndPoints)/2;
    hold on
end
xl = xlim;
xlim([xl(1) b(nPred).XEndPoints+0.1])
xticks(tks)
xticklabels(voc_names(vocFound))
xtickangle(45)
ax = gca;
set(gca,'yscale','log')
ax.FontSize = ft;
ax.FontWeight = 'bold';
ylabel('prevalence (log scale)')
yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = 0.05 + p(1); 
set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(e)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')

median(predictFreqVocInfer(~isinf(predictFreqVocInfer)))
median(predictPrevVocInfer(~isinf(predictPrevVocInfer)))
[c,p] = corr(numSeq',sum(~isinf(predictTime),2))
[c,p] = corr(numSeq',sum(predictTime>0,2))
[c,p] = corr(numSeq',sum(predictTime1>0,2))

% exportgraphics(gcf,['results' filesep 'figures' filesep 'inferSummary.png'],'Resolution',600)

