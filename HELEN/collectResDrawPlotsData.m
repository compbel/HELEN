%% 
% load basic data
clear;

continentsAll = {'Europe','Europe','North_America','North_America','Asia','Asia','Africa','Africa','Oceania','South_America','South_America','Europe','Europe','Europe','Europe','Europe'};
countriesAll = {'United_Kingdom','Germany','USA','Canada','Japan','India','South Africa','Kenya','Australia','Brazil','Peru','France','Sweden','Spain','Denmark','Italy'};
countriesAllPrint = countriesAll;
countriesAllPrint{1} = 'UK';
countriesAllPrint{7} = 'RSA';
% countriesAllPrint{10} = 'NZ';

% country = countriesAll{3}

% load(['data_' country '_2022.mat']);
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
% basic properties of data

numSeq = zeros(1,length(countriesAll));
for c = 1:length(countriesAll)
    c
    load(['data_' countriesAll{c} '_2022.mat']);
    indSeq = (sampTimeSeqD <= time_test_end);
    numSeq(c) = sum(indSeq);
end
(numSeq(1)+numSeq(3))/sum(numSeq)

lgComp = zeros(length(countriesAll),length(time_points));
lgComp2 = zeros(length(countriesAll),length(time_points));
bestTestDeg= zeros(length(countriesAll),length(time_points));

for c = 1:length(countriesAll)
%     c = 2
    country = countriesAll{c}
    degDistrFile = ['results' filesep 'network structure' filesep country '_aic_bic.mat'];
%     degCorrFile = ['results' filesep 'network structure' filesep country '_correlation.mat'];
    load(degDistrFile)
%     load(degCorrFile)
    
    for t = 1:length(time_points)
        outdir = ['alleles_joint_lm_2022_' country filesep char(time_points(t))];
        load([outdir filesep 'edges_1_1.mat']);        
        G = graph(E(:,1),E(:,2));
        if ~isempty(E)
            lgNode = max(max(E));
        else
            lgNode = 0;
        end
        if lgNode < nMut
            G = addnode(G,nMut - lgNode);
        end
        [comps,compsizes] = conncomp(G);
        [maxSizes,imaxSizes] = maxk(compsizes,2);
        lgComp(c,t) = maxSizes(1)/sum(compsizes);
        lgComp2(c,t) = maxSizes(2)/sum(compsizes);
    end
    

%     Dsize = readtable(netFileSize);
%     Ddeg = readtable(degFileSize);
%     T = find(~isnan(Dsize.largestComponentRelativeSize),1,'last');
%     figure
%     plot(Dsize.timepoints(1:T),Dsize.largestComponentRelativeSize(1:T),Dsize.timepoints(1:T),Dsize.secondLargestComponentRelativeSize(1:T),'LineWidth',2);
%     title('Sizes of the largest and second largest connected components')
%     ylabel('Relative size','fontweight','bold','fontsize',12)
%     legend('Largest component','Second largest component','fontweight','bold','fontsize',12,'Location','northwest')
%     a = get(gca,'XTickLabel');  
%     set(gca,'XTickLabel',a,'FontWeight','bold')
%     exportgraphics(gca,['results' filesep 'Uk_relative_size.png'],'Resolution',300)

    bic = AIC_BIC_mat(:,2:2:size(AIC_BIC_mat,2));
    [maxT,iMaxT] = min(bic,[],2);
    bestTestDeg(c,:) = iMaxT;
    
%     figure
%     plot(time_points,bic,'LineWidth',2)
%     grid on
%     title('(C) BIC for degree distribution fit')
%     ylabel('BIC','fontweight','bold','fontsize',12)
%     legend('Negative Binomial','Beta Negative Binomial','Poisson','Yule-Simon','Pareto','Generalized Pareto','fontweight','bold','fontsize',12,'Location','northwest')
%     a = get(gca,'XTickLabel');  
%     set(gca,'XTickLabel',a,'FontWeight','bold')
%     exportgraphics(gca,['results' filesep 'figures' filesep country '_BIC_deg.png'],'Resolution',300)

    % figure
    % title('(D) Correlation between network degree and VOC/VOI degree')
    % grid on
    % yyaxis left
    % plot(time_points,corrsP,'LineWidth',2);
    % ylabel('Correlation','fontweight','bold','fontsize',12)
    % ylim([0 1])
    % a = get(gca,'XTickLabel');  
    % set(gca,'XTickLabel',a,'FontWeight','bold')
    % yyaxis right
    % plot(time_points,pvaluesP,'LineWidth',2);
    % ylabel('p-value','fontweight','bold','fontsize',12)
    % ylim([-0.1 1])
    % ylim('tickaligned')
    % exportgraphics(gca,['results' filesep 'figures' filesep country '_corr_deg.png'],'Resolution',300)
end

[h,p] = kstest2(lgComp(:),lgComp2(:))
histcounts(bestTestDeg(:))/numel(bestTestDeg)
median(arrayfun(@(x)find(bestTestDeg(x,:)~=4,1,'last'),1:size(bestTestDeg,1)))
%% 

fig = figure;
tile = tiledlayout(2,3);
set(gcf, 'Position',  [100, 100, 1600, 1000])

[~,ind] = sort(numSeq,'descend');
names = categorical(countriesAllPrint(ind));
names = reordercats(names,countriesAllPrint(ind));
nexttile
bar(names,numSeq(ind))
% title('(A) Dataset sizes for analyzed countries')
grid on
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'FontWeight','bold')
ylabel('number of sequences','fontweight','bold')
xl = xlim;
yl = ylim;
text(length(names),0.99*yl(2),'(a)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')

c=3; 
t=18;
v = 5;
time_points(t)
outdir = ['alleles_joint_lm_2022_' countriesAll{c} filesep char(time_points(t))];
load([outdir filesep 'edges_1_1.mat']);        
G = graph(E(:,1),E(:,2));
lgNode = max(max(E));
if lgNode < nMut
    G = addnode(G,nMut - lgNode);
end
[comps,compsizes] = conncomp(G);
ind = (comps == find(compsizes==max(compsizes)));
H = subgraph(G,ind);
voc_comp = voc(:,ind);
nexttile([2 2]);
h = plot(H,'layout','force','EdgeColor','k');
highlight(h,find(voc_comp(v,:)),'NodeColor',voc_colors(v,:),'MarkerSize',4);
xl = xlim;
yl = ylim;
text(0.79*xl(2),0.79*yl(2),'(c)','FontSize',20,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
hold on
voc_leg = plot(NaN,NaN,'Color',voc_colors(v,:),'LineWidth',4);
hold off
legend(voc_leg,voc_names{v},'location', 'Southwest','FontSize',12,'FontWeight','bold')

nexttile
plot(time_points,median(lgComp,1),'b',time_points,min(lgComp,[],1),'b--',time_points,max(lgComp,[],1),'b--','LineWidth',2);
grid on
hold on
plot(time_points,median(lgComp2,1),'r',time_points,min(lgComp2,[],1),'r--',time_points,max(lgComp2,[],1),'r--','LineWidth',2);
grid on
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'FontWeight','bold')
ylabel('component size, %','fontweight','bold')
xl = xlim;
yl = ylim;
tile.TileSpacing = 'tight';
text(max(time_points)+days(60),yl(2),'(b)','FontSize',18,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
hold on
voc_leg(1) = plot(NaN,NaN,'Color','b','LineWidth',4);
voc_leg(2) = plot(NaN,NaN,'Color','r','LineWidth',4);
hold off
legend(voc_leg,{'Largest component','Second largest component'},'location', 'Northwest','FontSize',12)


exportgraphics(fig,['results' filesep 'figures' filesep 'data_net1.png'],'Resolution',600)


%% 
% save voc-related data
M_all = [];
sampTimeSeqD_all = [];
countrySeq_all = [];
for c = 1:length(countriesAll)
    c
    freq_voc_all = zeros(size(voc,1),length(time_points));
    cumInc_voc_all = zeros(size(voc,1),length(time_points));
    relPrev_voc_all = zeros(size(voc,1),length(time_points));
    
    load(['data_' countriesAll{c} '_2022.mat']);
    outdirRes = ['results' filesep 'pvalues' filesep countriesAll{c} '_pvalues'];
    for t = 1:length(time_points)
        [c t]
        for s = 1:size(voc,1)

            outfile = [outdirRes filesep 'res_s' int2str(s) '_t' int2str(t) '.mat'];
            if ~exist(outfile,'file')
                continue;
            end
            load(outfile)
            cumInc_voc_all(s,t) = cumInc_voc;
            freq_voc_all(s,t) = freq_voc;
            ind = (sum(M(:,voc(s,:)),2) == sum(voc(s,:)))';
            if t >= 2    
                ind1 = (sampTimeSeqD > time_points(t-1)) & (sampTimeSeqD <= time_points(t));
                relPrev_voc_all(s,t) = sum(ind&ind1)/sum(ind1);
            else
                relPrev_voc_all(s,t) = freq_voc_all(s,t);
            end
        end
    end
    save(['vocs_' countriesAll{c} '.mat'],'cumInc_voc_all','freq_voc_all','relPrev_voc_all');
    M_all = [M_all; M];
    sampTimeSeqD_all = [sampTimeSeqD_all sampTimeSeqD];
    countrySeq_all = [countrySeq_all c*ones(1,length(sampTimeSeqD))];
end
save(['data_all_2022.mat'],'M_all','sampTimeSeqD_all','countrySeq_all','-v7.3');

%% 
% add new stats to inferred voc data
for c = 1:length(countriesAll)
    res_file = ['results' filesep countriesAll{c} '_densestFixSize.mat'];
    load(res_file);
    load(['data_' countriesAll{c} '_2022.mat']);
    for t = 1:length(time_points)
            [c t]
            vocInferCummFreq = zeros(1,length(stats{t}.vocInfer));
            vocInferPrev = zeros(1,length(stats{t}.vocInfer));
            for l = 1:length(stats{t}.vocInfer)
                ind1 = (sum(M(:,stats{t}.vocInfer{l}),2) == length(stats{t}.vocInfer{l}))';
                ind2 = (sampTimeSeqD <= time_points(t));
                vocInferCummFreq(l) = sum(ind1&ind2)/sum(ind2);
                if t >=2
                    ind3 = (sampTimeSeqD > time_points(t-1));
                    vocInferPrev(l) = sum(ind1&ind2&ind3)/sum(ind2&ind3);
                end
            end
            stats{t}.vocInferCummFreq = vocInferCummFreq;
            stats{t}.vocInferPrev = vocInferPrev;
    end
    save(['results' filesep country '_densestFixSizeExt.mat'],'stats');
end
