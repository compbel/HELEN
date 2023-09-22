%% 
% load basic data
clear;

continentsAll = {'Europe','Europe','North_America','North_America','Asia','Asia','Africa','Africa','Oceania','South_America','South_America','Europe','Europe','Europe','Europe','Europe'};
countriesAll = {'United_Kingdom','Germany','USA','Canada','Japan','India','South_Africa','Kenya','Australia','Brazil','Peru','France','Sweden','Spain','Denmark','Italy'};
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
voc_design_fig{4} = datetime('05-11-2021','InputFormat','MM-dd-yyyy');
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
% time_test_end = datetime('11-01-2021','InputFormat','MM-dd-yyyy');
time_test_end = datetime('12-31-2021','InputFormat','MM-dd-yyyy');

time_test_end_extend = datetime('04-22-2022','InputFormat','MM-dd-yyyy');

time_points = time_test_begin:caldays(15):time_test_end;
time_points_extend = time_test_begin:caldays(15):time_test_end_extend;

dataset = '1';
outfolder = ['HELEN_data' dataset];
voc = voc(1:10,:);

%% 
% densest subgraph analysis
moveLegend = zeros(1,length(countriesAll));
% moveLegend([3 4 7 10 11]) = 1;
moveLegend([10]) = 1;
predictLevels = [0.8 0.9];
predictInterv = cell(1,length(predictLevels));
predictInterv1 = cell(1,length(predictLevels));
predictTime = cell(1,length(predictLevels)); 
predictMetr = cell(1,length(predictLevels));
predictVocFreq = cell(1,length(predictLevels));
predictVocPrev = cell(1,length(predictLevels));
allF = zeros(length(countriesAll),length(time_points));
for l = 1:length(predictLevels)
    predictInterv{l} = -Inf*ones(length(countriesAll),size(voc,1));
    predictInterv1{l} = -Inf*ones(length(countriesAll),size(voc,1));
    predictTime{l} = NaT(length(countriesAll),size(voc,1));
    predictMetr{l} = zeros(length(countriesAll),size(voc,1));
    predictVocFreq{l} = -Inf*ones(length(countriesAll),size(voc,1));
    predictVocPrev{l} = -Inf*ones(length(countriesAll),size(voc,1));
end
scatterPointsF = [];
scatterPointsPrevDensest = [];
scatterPointsPrev = [];
scatterPointsFreq = [];
scatterPointsT = repmat(time_points,1,length(countriesAll));
scatterPointsV = [];
scatterPoints1Perc = [];

% bar plot f-scores per country
nRow = 4;
nCol = 4;
figure
[ha, pos] = tight_subplot(nRow,nCol,[.035 .025],[.1 0.05],[.1 .1]);
set(gcf, 'Position',  [100, 100, 2650, 1250])
% set(gcf, 'Position',  [100, 100, 1250 2000])

count = 0;
for c = 1:length(countriesAll)
    count = count + 1;
    axes(ha(count)); 
    countriesAll{c}
    subgr_file = [countriesAll{c} '_densest_subgraph.mat'];
    load([outfolder filesep 'Densest_subgraphs' filesep subgr_file]);
    % load([outfolder filesep 'Preprocessed_input_data' filesep 'vocs_' countriesAll{c} '.mat']);
    load([outfolder filesep 'Preprocessed_input_data' filesep 'vocs_' countriesAll{c} '_GIS.mat']);

    
    
    [prc,tm1perc] = max(relPrev_voc_all>=0.01,[],2);
    tm1perc(prc == 0) = length(time_points);
    
    metric = f_densest;
    [best_metric,best_voc_tp] = max(metric(:,1:10),[],2);
    rel_inc_best_voc_tp = relPrev_voc_all(sub2ind(size(relPrev_voc_all),best_voc_tp',1:length(time_points)));
    freq_best_voc_tp = freq_voc_all(sub2ind(size(freq_voc_all),best_voc_tp',1:length(time_points)));
    
    allF(c,:) = best_metric';
    scatterPointsF = [scatterPointsF best_metric'];
    scatterPointsPrev = [scatterPointsPrev rel_inc_best_voc_tp];
    % scatterPointsPrevDensest = [scatterPointsPrevDensest prev_densest'];
    scatterPointsFreq = [scatterPointsFreq freq_best_voc_tp];
    scatterPointsV = [scatterPointsV best_voc_tp'];
    scatterPoints1Perc = [scatterPoints1Perc time_points(tm1perc(best_voc_tp))];
    
    best_voc_freq = zeros(1,length(time_points));
    best_voc_prev = zeros(1,length(time_points));
    for i = 1:length(time_points)
        best_voc_prev(i) = relPrev_voc_all(best_voc_tp(i),i);
        best_voc_freq(i) = freq_voc_all(best_voc_tp(i),i);
    end
%     predictVocFreq(c,:) = best_voc_freq;
    
    vocPredict = (unique(best_voc_tp))';
    for v = vocPredict
        for l = 1:length(predictLevels)
            t = find((best_voc_tp == v)&(round(best_metric,2) >= predictLevels(l)),1,'first');
            if ~isempty(t)
                predictInterv{l}(c,v) = -days(time_points(t) - voc_design{v});
                predictTime{l}(c,v) = time_points(t);
                predictMetr{l}(c,v) = best_metric(t);
                predictVocFreq{l}(c,v) = best_voc_freq(t);
                predictVocPrev{l}(c,v) = best_voc_prev(t);

                % tm1perc = time_points(min([find(relPrev_voc_all(v,:)>=0.01,1,'first') length(time_points)]));
                % if tm1perc == time_points(end)
                %     tm1perc = max([tm1perc voc_design{v}]);
                % end
                tm1perc = find(relPrev_voc_all(v,:)>=0.01,1,'first');
                if isempty(tm1perc)
                    tm1perc = find(relPrev_voc_all(v,:) == max(relPrev_voc_all(v,:)),1,'first');
                end
                predictInterv1{l}(c,v) = -days(time_points(t)-time_points(tm1perc));
            end
        end
    end
    
    % figure
    b = bar(time_points,best_metric,0.5,'FaceColor','flat');
    
    title(strrep(countriesAllPrint{c},'_',' '))
    % xlim([time_test_begin-days(15) voc_design{5}+days(15)])
    grid on
    for i = 1:length(time_points)
        b.CData(i,:) = voc_colors(best_voc_tp(i),:);
    end
    ylabel('voc-specific f-score','fontweight','bold','fontsize',12)
    ylim([0 1])
%     if max(ytips) >= 0.97
%         ylim([0 1.1])
%     end
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'FontWeight','bold')
    vocs_densest = (unique(best_voc_tp))';
    for s = vocs_densest
        xl = xline(voc_design_fig{s},'--','Color',voc_colors(s,:),'LineWidth',2);
    %     xl.LabelOrientation = 'horizontal';
    end
    hold on
    voc_leg = [];
    for i = 1:length(vocs_densest)
        voc_leg(i) = plot(NaN,NaN,'Color',voc_colors(vocs_densest(i),:),'LineWidth',2);
    end
    hold off
    if moveLegend(c)
        legend(voc_leg,voc_names(vocs_densest),'location', 'northeast')
    else
        legend(voc_leg,voc_names(vocs_densest),'location', 'northwest')
    end
    % exportgraphics(gca,['Figures' filesep countriesAll{c} '_densestFscore_D' dataset '.png'],'Resolution',300)
end
exportgraphics(gcf,['Figures' filesep 'densestFscore_D' dataset '.png'],'Resolution',600)

%% 
% bar plot f-scores per voc
ft = 11;
nVoc = 5;
fig = figure
tile = tiledlayout(nVoc,1);
axes = [];
for v = 1:nVoc
    data = zeros(length(time_points),length(countriesAll));
    toFill = zeros(length(time_points),length(countriesAll));
    for t = 1:length(time_points)
        ind = find(scatterPointsT==time_points(t)&(scatterPointsV==v));
        [data(t,1:length(ind)),indSort] = sort(scatterPointsF(ind),'descend');
%         toFill(t,1:length(ind)) = scatterPoints1Perc(ind(indSort)) >= time_points(t);
        % toFill(t,1:length(ind)) = scatterPointsPrevDensest(ind(indSort)) <= 0.01;
    end
    k = find(sum(data,1)>0,1,'last');
    data = data(:,1:k);
    toFill = toFill(:,1:k);
%     figure
    axes(v) = nexttile;
%     b = bar(axes(v),time_points,data,'BarWidth', 1,'FaceColor','flat','EdgeColor',voc_colors(v,:));
    b = bar(axes(v),time_points,data,'BarWidth', 1,'FaceColor','flat','EdgeColor','k','LineWidth',0.25);
    grid on
    for i = 1:length(b)
        for j = 1:size(data,1)
            b(i).CData(j,:) = voc_colors(v,:);
%             if toFill(j,i) == 1
%                 b(i).CData(j,:) = voc_colors(v,:);
%             else
%                 b(i).CData(j,:) = [1 1 1];
%             end
        end
    end
%     ylabel('f-score','fontweight','bold','fontsize',14)
%     ylim([0 1])
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'FontWeight','bold','FontSize',ft)
    set(gcf, 'Position',  [100, 100, 1900, 650])
    xl = xline(voc_design_fig{v},'--','Color','k','LineWidth',3);
    yl = yline(0.8,':','Color','k','LineWidth',0.5);
    if v ~= nVoc
        xticklabels(axes(v),{});
    end
    if v == 1
        hold on
        voc_leg = [];
        for i = 1:nVoc
            voc_leg(i) = plot(NaN,NaN,'Color',voc_colors(i,:),'LineWidth',4);
        end
        hold off
        legend(voc_leg,voc_names,'location', 'Northeast')
    end
    yticks([0 0.5 0.8]);
end
linkaxes(axes,'x');
% title(tile,'My Title')
% xlabel(tile,'time')
ylabel(tile,'fscore','fontweight','bold','fontsize',ft)
xl = xlim;
% xlim([xl(1) voc_design{5} + days(7)]) 
tile.TileSpacing = 'none';
exportgraphics(fig,['Figures' filesep 'densestFscoreAll_D' dataset '.png'],'Resolution',600)

%% 


totalDetect = zeros(1,5);
prevDetect = [];
freqDetect = [];
for s = 1:5
    ind = (scatterPointsV == s) & (scatterPointsT <= voc_design{s}) & (scatterPointsF >= 0.8);
    totalDetect(s) = sum(ind);
    ind = find((scatterPointsV == s) & (scatterPointsF >= 0.8));
    ind_tp = ind(scatterPointsT(ind)==min(scatterPointsT(ind)));
    prevDetect = [prevDetect scatterPointsPrev(ind_tp)];
    freqDetect = [freqDetect scatterPointsFreq(ind_tp)];
end
mean(prevDetect)
median(prevDetect)
std(prevDetect)
mean(freqDetect)
median(freqDetect)
std(freqDetect)

sum(scatterPointsF >= 0.8)/length(scatterPointsF)
sum(totalDetect)/sum(scatterPointsF >= 0.8)
predictTime = predictInterv{1};
predictTime1 = predictInterv1{1};
predictVocFreq = predictVocFreq{1};
predictVocPrev = predictVocPrev{1};

%% 
% inferred communities forecasting depths and frequency bar plot

ft = 11;
vocFound = find(sum(~isinf(predictTime),1)>0);
w0 = 4;

nRow = 1;
nCol = 4;
figure
[ha, pos] = tight_subplot(nRow,nCol,[.045 .03],[.1 0.05],[.1 .1]);
set(gcf, 'Position',  [100, 100, 2900, 770])

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
% yh = get(gca,'ylabel');
% p = get(yh,'position');
% p(1) = 0.05 + p(1); 
% set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(a)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')

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
% yh = get(gca,'ylabel');
% p = get(yh,'position');
% p(1) = 0.2 + p(1); 
% set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(b)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')


% figure
axes(ha(3));
title('(c) Cummulative frequencies at times of variant calling')
% set(gcf, 'Position',  [100, 100, 1190, 805])
xNext = 1;
tks = zeros(1,length(vocFound));
for v = 1:length(vocFound)
    ptSorted = sort(predictVocFreq(:,vocFound(v)),'descend');
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
% yh = get(gca,'ylabel');
% p = get(yh,'position');
% p(1) = 0.05 + p(1); 
% set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(c)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')

% figure
axes(ha(4));
title('(d) Prevalences at times of variant calling')
% set(gcf, 'Position',  [100, 100, 1190, 805])
xNext = 1;
tks = zeros(1,length(vocFound));
for v = 1:length(vocFound)
    ptSorted = sort(predictVocPrev(:,vocFound(v)),'descend');
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
% yh = get(gca,'ylabel');
% p = get(yh,'position');
% p(1) = 0.05 + p(1); 
% set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(d)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')

exportgraphics(gcf,['Figures' filesep 'densestSummary_D' dataset '.png'],'Resolution',600)

%% 

vocFound = find(sum(~isinf(predictTime),1)>0)
sum(~isinf(predictTime),1)

sum(predictTime >= 0,1)
sum(predictTime1 >= 0,1)

sum(scatterPointsF >= 0.8)/length(scatterPointsF)
sum(~isinf(predictTime(:)))/numel(predictTime)
sum(sum(~isinf(predictTime(:,1:5))))/numel(predictTime(:,1:5))

median(predictVocFreq(~isinf(predictVocFreq)))
median(predictVocPrev(~isinf(predictVocPrev)))

median(predictVocFreq(predictTime1>0 | predictTime>0))
median(predictVocPrev(predictTime1>0 | predictTime>0))

sum(predictTime(:)>0)/sum(~isinf(predictTime(:)))
sum(predictTime1(:)>0)/sum(~isinf(predictTime1(:)))

median(predictTime(~isinf(predictTime)))
median(predictTime1(~isinf(predictTime1)))

median(predictTime(predictTime>0))
median(predictTime1(predictTime1>0))

max(predictTime,[],1)
max(predictTime1,[],1)