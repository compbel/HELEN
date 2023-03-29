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
voc_design_fig{2} = datetime('12-18-2020','InputFormat','MM-dd-yyyy')+days(5);
voc_design_fig{3} = datetime('01-11-2021','InputFormat','MM-dd-yyyy')+days(6);
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
time_test_end = datetime('11-01-2021','InputFormat','MM-dd-yyyy');
time_test_end_extend = datetime('04-22-2022','InputFormat','MM-dd-yyyy');
time_points = time_test_begin:caldays(15):time_test_end;
time_points_extend = time_test_begin:caldays(15):time_test_end_extend;

%% 
% samplng-based p-values figures separately

lw = 4;
correl = zeros(length(countriesAll),size(voc,1));
pvalue = zeros(length(countriesAll),size(voc,1));
hvalue = zeros(length(countriesAll),size(voc,1));
correl1 = zeros(length(countriesAll),size(voc,1));
pvalue1 = zeros(length(countriesAll),size(voc,1));
granger = zeros(length(countriesAll),size(voc,1));
optLags = zeros(length(countriesAll),size(voc,1));
optLags1 = zeros(length(countriesAll),size(voc,1));
pgranger = zeros(length(countriesAll),size(voc,1));
predictTime = -Inf*ones(length(countriesAll),size(voc,1));
predictTime1 = -Inf*ones(length(countriesAll),size(voc,1));
times1Perc = NaT(length(countriesAll),size(voc,1));
predictInc = -Inf*ones(length(countriesAll),size(voc,1));
predictFreq = -Inf*ones(length(countriesAll),size(voc,1));
predictPrev = -Inf*ones(length(countriesAll),size(voc,1)-3);
predictValue = -Inf*ones(length(countriesAll),size(voc,1));
padj_voc_all = cell(1,length(countriesAll));

for c = 1:length(countriesAll)
    countriesAll{c}
    p_voc_time = zeros(size(voc,1),length(time_points));
    padj_voc_all{c} = zeros(size(voc,1),length(time_points));
    perc_lgcomp_time = zeros(size(voc,1),length(time_points));
    outdirRes = ['HELEN data' filesep 'p-values' filesep countriesAll{c}];
    load(['HELEN data' filesep 'Preprocessed input data' filesep 'vocs_' countriesAll{c} '.mat']);

    for t = 1:length(time_points)
        for s = 1:size(voc,1)

            outfile = [outdirRes filesep 'res_s' int2str(s) '_t' int2str(t) '.mat'];
            load(outfile)
            p_voc_time(s,t) = p_voc;
            perc_lgcomp_time(s,t) = perc_lgcomp_voc;
        end
    end
    relPrev_voc_all(isnan(relPrev_voc_all)) = 0;
    
    for s = 1:(size(voc,1))
        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_voc_time(s,:));
        ptime = find((h==1)&(perc_lgcomp_time(s,:)>=0.8),1,'first');
        padj_voc_all{c}(s,:) = adj_p;
        times1Perc(c,s) = time_points(min([find(relPrev_voc_all(s,:)>=0.01,1,'first') length(time_points)]));
         
        

        if ~isempty(ptime)
            predictTime(c,s) = days(voc_design{s} - time_points(ptime));
            predictInc(c,s) = relPrev_voc_all(s,ptime);
            predictFreq(c,s) = freq_voc_all(s,ptime);
            predictPrev(c,s) = relPrev_voc_all(s,ptime);
            predictValue(c,s) = p_voc_time(s,ptime);
            
            tm1perc = min([find(relPrev_voc_all(s,:)>=0.01,1,'first') length(time_points)]);
            predictTime1(c,s) = days(time_points(tm1perc)- time_points(ptime));
            times1Perc(c,s) = time_points(tm1perc);
        end

        if sum(relPrev_voc_all(s,:)) == 0
            relPrev_voc_all(s,:) = linspace(0.0001,0.001,length(time_points));
        end
        [~,mi] = max(relPrev_voc_all(s,:));
        [correl1(c,s),pvalue1(c,s)] = corr(padj_voc_all{c}(s,1:mi)',relPrev_voc_all(s,1:mi)','Type','Pearson');

        [xcf,lags,bounds] = crosscorr(padj_voc_all{c}(s,1:mi)',relPrev_voc_all(s,1:mi)');
        ind = (xcf >= bounds(1)) | (xcf <= bounds(2));
        xcf = xcf(ind);
        lags = lags(ind);
        if sum(~isnan(xcf))>0
            correl(c,s) = xcf(find(abs(xcf)==max(abs(xcf)),1,'last'));
            optLags(c,s) = lags(find(abs(xcf)==max(abs(xcf)),1,'last'));
            hvalue(c,s) = 1;
        else
            correl(c,s) = NaN;
            optLags(c,s) = NaN;
            pvalue(c,s) = NaN;
        end
    end
end

k = 5;
correl = correl(:,1:k);
pvalue = pvalue(:,1:k);
hvalue = hvalue(:,1:k);
optLags = optLags(:,1:k);
ind = hvalue(:)==1;
sum((hvalue(:)==1)&(correl(:)<=-0.2))/numel(correl)
signCorrel = sort(correl(ind));
signLags = sort(optLags(ind));
signCorrel = signCorrel(ceil(0.05*length(signCorrel)):(end - ceil(0.05*length(signCorrel))+1));
signLags = signLags(ceil(0.05*length(signLags)):(end - ceil(0.05*length(signLags))+1));
median(signCorrel)
median(signLags)
[signCorrel(1) signCorrel(end)]
[signLags(1) signLags(end)]


scatterPointsY = [];
scatterPointsY1 = [];
scatterPointsY2 = [];
scatterPointsY3 = [];
scatterPointsX = [];
indMax = [];
indMax1 = [];
for s = 1:(size(voc,1)-3)
    scatterPointsX = [scatterPointsX s*ones(1,sum(~isinf(predictTime(:,s))))];
    scatterPointsY1 = [scatterPointsY1 (predictInc(~isinf(predictInc(:,s)),s))'];
    scatterPointsY2 = [scatterPointsY2 (predictFreq(~isinf(predictFreq(:,s)),s))'];
    tms = (unique(predictTime(~isinf(predictTime(:,s)),s)))';
    for t = tms
        ind = find(predictTime(:,s) == t);
        for j = 1:length(ind)
            scatterPointsY = [scatterPointsY predictTime(ind(j),s)-10*(j-1)];
        end
    end
    tms1 = (unique(predictTime1(~isinf(predictTime1(:,s)),s)))';
    for t = tms1
        ind = find(predictTime1(:,s) == t);
        for j = 1:length(ind)
            scatterPointsY3 = [scatterPointsY3 predictTime1(ind(j),s)-10*(j-1)];
        end
    end
    try
        indMax = [indMax find((scatterPointsX == s)&(scatterPointsY == max(scatterPointsY(scatterPointsX == s))))];
        indMax1 = [indMax1 find((scatterPointsX == s)&(scatterPointsY3 == max(scatterPointsY3(scatterPointsX == s))))];
    catch
        [];
    end
end

%% 
% samplng-based p-values forecasting depths and frequency bar plot
predictTime = predictTime(:,1:10);
predictTime1 = predictTime1(:,1:10);
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
% title('(b) Forecasting depth wrt 1% prevalence time')
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
ylabel('FD^{prev}, days')
yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = 0.25 + p(1); 
set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(b)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
% exportgraphics(gcf,['results' filesep 'figures' filesep 'p_predict_freq.png'],'Resolution',600)

% figure
axes(ha(2)); 
title('(c) Forecasting depth wrt designation time')
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
ylabel('FD^{des}, days')
yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = 0.2 + p(1); 
set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(c)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
% exportgraphics(gcf,['results' filesep 'figures' filesep 'p_predict_des.png'],'Resolution',600)

% exportgraphics(gcf,['results' filesep 'figures' filesep 'p_freq.png'],'Resolution',300)
% figure
axes(ha(3)); 
title('(d) Cummulative frequencies at times of variant calling')
% set(gcf, 'Position',  [100, 100, 1190, 805])
xNext = 1;
tks = zeros(1,length(vocFound));
for v = 1:length(vocFound)
%     ptSorted = sort(predictTime(~isinf(predictTime(:,v)),v),'descend');
    ptSorted = sort(predictFreq(:,vocFound(v)),'descend');
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
ylabel('cummulative frequency (log scale)');
yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = 0.05 + p(1); 
set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(d)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
% exportgraphics(gcf,['results' filesep 'figures' filesep 'p_freq.png'],'Resolution',300)

% figure
axes(ha(4)); 
title('(e) Prevalences at times of variant calling')
% set(gcf, 'Position',  [100, 100, 1190, 805])
xNext = 1;
tks = zeros(1,length(vocFound));
for v = 1:length(vocFound)
    ptSorted = sort(predictPrev(:,vocFound(v)),'descend');
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
ylabel('prevalence (log scale)');
yh = get(gca,'ylabel');
p = get(yh,'position');
p(1) = 0.05 + p(1); 
set(yh,'position',p);
xl = xlim;
yl = ylim;
text(0.99*xl(2),0.99*yl(2),'(e)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')


exportgraphics(gcf,['results' filesep 'figures' filesep 'pvalueSummary1.png'],'Resolution',600)

sum(sum(~isinf(predictTime1)))/numel(predictTime1)
sum(sum(~isinf(predictTime1(:,1:5))))/numel(predictTime1(:,1:5))
sum(sum(predictTime>0))/sum(sum(~isinf(predictTime)))
sum(sum(predictTime1>0))/sum(sum(~isinf(predictTime1)))


predictFreq = predictFreq(:,1:10);
aux = predictFreq(~isinf(predictFreq));
median(aux)

predictPrev = predictPrev(:,1:10);
aux = predictPrev(~isinf(predictPrev));
median(aux)

aux = predictTime(predictTime>0);
median(aux)

aux = predictTime1(predictTime1>0);
median(aux)

[c,p] = corr(numSeq',sum(~isinf(predictTime(:,1:10)),2))
%% 

% samplng-based p-values figures together
lw = 2;
ft = 8;
p_voc_all = cell(1,length(countriesAll));
perc_lgcomp_all = cell(1,length(countriesAll));
relInc_all = cell(1,length(countriesAll));
for c = 1:length(countriesAll)
    countriesAll{c}
    p_voc_time = zeros(size(voc,1),length(time_points));
    perc_lgcomp_time = zeros(size(voc,1),length(time_points));
    
    outdirRes = ['HELEN data' filesep 'p-values' filesep countriesAll{c}];
    load(['HELEN data' filesep 'Preprocessed input data' filesep 'vocs_' countriesAll{c} '.mat']);

    for t = 1:length(time_points)
        for s = 1:size(voc,1)

            outfile = [outdirRes filesep 'res_s' int2str(s) '_t' int2str(t) '.mat'];
            if ~exist(outfile,'file')
                continue;
            end
            load(outfile)
            p_voc_time(s,t) = p_voc;
            perc_lgcomp_time(s,t) = perc_lgcomp_voc;
        end
    end
    p_voc_all{c} = p_voc_time;
    perc_lgcomp_all{c} = perc_lgcomp_time;
    relInc_all{c} = relPrev_voc_all;
end 
%% 

% vocs per country
lw = 3.5;
ft = 11;
for c = 1:length(countriesAll) 
    nRow = 2;
    nCol = 4;
    [ha, pos] = tight_subplot(nRow,nCol,[.045 .03],[.1 0.05],[.1 .1]);
    %     [~,excl] = mink(cumInc_voc_all(6:end,end),4);
    [~,excl] = maxk(padj_voc_all{c}(6:end,end),size(voc,1) - nRow*nCol);
    excl = excl + 5;

    set(gcf, 'Position',  [100, 100, 1400, 600])
    count = 1;
    for s = 1:(size(voc,1))
        if ~isempty(find(excl == s))
            continue;
        end
        axes(ha(count)); 
        title([voc_names{s} ' ' countriesAllPrint{c}],'fontsize',ft)
        grid on
        yyaxis left
        plot(time_points,padj_voc_all{c}(s,:),'LineWidth',lw);
        ylim([0 1])
        ax=gca;
        ax.FontSize = ft;
        if mod(count,nCol) == 1
            ylabel('p-value','fontweight','bold','fontsize',ft)
        else
            set(gca,'YTickLabel',[]);
        end
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'FontWeight','bold')
        if count <= (nRow-1)*nCol
            set(gca,'XTickLabel',[]);
        end
        yyaxis right
        plot(time_points,relInc_all{c}(s,:),'LineWidth',lw);
        if mod(count,nCol) == 0
            ylabel('prevalence','fontweight','bold','fontsize',ft)
        end
        if max(relInc_all{c}(s,:) > 0)
            ylim('tickaligned')
        else
            ylim([-0.1 1])
    %         ylim('tickaligned')
        end
        if ~isempty(voc_design{s})
            xline(voc_design_fig{s},'--k','LineWidth',lw);
        end
%         ptime = find(p_voc_all(s,:)<=0.05,1,'first');
        ptime = find((padj_voc_all{c}(s,:)<=0.05)&(perc_lgcomp_all{c}(s,:)>=0.8),1,'first');
        if ~isempty(ptime)
            xline(time_points(ptime),'--m','LineWidth',lw);
        end
        xline(times1Perc(c,s),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',lw);
        xl = xlim;
        yl = ylim;
        if count == 1
            text(time_points(8),yl(2),'(a)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
        end
        count = count+1;
    end
%     exportgraphics(gcf,['results' filesep 'figures' filesep countriesAll{c} '_pvalues.png'],'Resolution',600)
end 
%% 

% countries per voc
lw = 2;
ft = 8;
for s = 1:(size(voc,1))  
    nRow = 4;
    nCol = 4;
    figure
    [ha, pos] = tight_subplot(nRow,nCol,[.035 .025],[.1 0.05],[.1 .1]);

    set(gcf, 'Position',  [100, 100, 1050, 900])
    count = 1;
    for c = 1:length(countriesAll) 
        axes(ha(count)); 
        title([voc_names{s} ' ' countriesAllPrint{c}],'fontsize',ft)
        grid on
        yyaxis left
        plot(time_points,padj_voc_all{c}(s,:),'LineWidth',lw);
        ylim([0 1])
        ax=gca;
        ax.FontSize = ft;
        if mod(count,nCol) == 1
            ylabel('p-value','fontweight','bold','fontsize',ft)
        else
            set(gca,'YTickLabel',[]);
        end
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'FontWeight','bold')
        if count <= (nRow-1)*nCol
            set(gca,'XTickLabel',[]);
        end
        yyaxis right
        plot(time_points,relInc_all{c}(s,:),'LineWidth',lw);
        if mod(count,nCol) == 0
            ylabel('prevalence','fontweight','bold','fontsize',ft)
        end
        if max(relInc_all{c}(s,:) > 0)
            ylim('tickaligned')
        else
            ylim([-0.1 1])
    %         ylim('tickaligned')
        end
        if ~isempty(voc_design{s})
            xline(voc_design{s},'--k','LineWidth',lw);
        end
%         ptime = find(p_voc_all(s,:)<=0.05,1,'first');
        ptime = find((padj_voc_all{c}(s,:)<=0.05)&(perc_lgcomp_all{c}(s,:)>=0.8),1,'first');
        if ~isempty(ptime)
            xline(time_points(ptime),'--m','LineWidth',lw);
        end
        xline(times1Perc(c,s),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',lw);
        count = count+1;
    end
%     exportgraphics(gcf,['results' filesep 'figures' filesep voc_names{s} '_pvalues_all.png'],'Resolution',600)
end

