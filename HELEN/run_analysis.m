clear;

%path to GUROBI LP solver
addpath 'C:\gurobi811\win64\matlab'
% addpath '/opt/gurobi810/linux64/matlab'


continents = {'Europe','Europe','North_America','North_America','Asia','Asia','Africa','Africa','Oceania','South_America','South_America','Europe','Europe','Europe','Europe','Europe'};
countries = {'United_Kingdom','Germany','USA','Canada','Japan','India','South Africa','Kenya','Australia','Brazil','Peru','France','Sweden','Spain','Denmark','Italy'};

time_test_begin = datetime('05-01-2020','InputFormat','MM-dd-yyyy');
time_test_end = datetime('11-01-2021','InputFormat','MM-dd-yyyy');
time_points = time_test_begin:caldays(15):time_test_end;

% voc/voi haplotypes
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
voc = false(10,nMut);
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

% days of first reported COVID-19 cases
day1 = {'30-Jan-2020','26-Jan-2020','21-Jan-2020','16-Feb-2020','21-Jan-2020','29-Jan-2020','04-Mar-2020','12-Mar-2020',...
    '23-Feb-2020','25-Feb-2020','05-Mar-2020','23-Jan-2020','31-Jan-2020','31-Jan-2020','26-Feb-2020','30-Jan-2020'};
day1 = datetime(day1,'InputFormat','dd-MMM-yyyy');

% Set the country to be analyzed
cnt = 7;
continent = continents{cnt};
country = countries{cnt};
day1 = day1(cnt);

%% 
% This section loads the data, prepares several data structures used by
% HELEN scripts and saves them for further use

% Load sequencing data and metadata downloaded from GISAID
seqF = fastaread(['Genomic data' filesep continent filesep 'sequences_' country '_protein.fasta']);
seq = char(seqF.Sequence);
seq = seq(:,1:end-1);

% Parse sampling times of sequences from GISAID data
sampTimeSeqD = NaT(1,length(seqF));
sampTimeSeq = zeros(1,length(seqF));
badDate = false(1,length(seqF));
formatDate = 'yyyy-MM-dd';
for i = 1:length(seqF)
    i
    C = strsplit(seqF(i).Header,'|');
    dateSeq = C{3};
    try
        z = datetime(dateSeq,'InputFormat',formatDate);
        sampTimeSeq(i) = days(z-day1);
        sampTimeSeqD(i) = z;
    catch
        badDate(i) = 1;
    end
end
ind = ~badDate;
seqF = seqF(ind);
seq = seq(ind,:);
sampTimeSeq = sampTimeSeq(ind);
sampTimeSeqD = sampTimeSeqD(ind);

% construct mutation matrix
refF = fastaread(['Genomic data' filesep 'reference_s_aa.fas']);
ref = char(refF.Sequence);
seq = fillGaps(seq,ref);
M = getMutMatr(seq,ref);
M = sparse(M);

% clear memory and save 
clear seq seqF aux badDate C dataEpid dateSeq ind sampCounts 
outdir = ['HELEN data1' filesep 'Preprocessed input data'];
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
save([outdir filesep 'data_' country '_2022.mat'],'-v7.3');

%%
% This section constructs epistatic networks and saves them for further use
load(['HELEN data1' filesep 'Preprocessed input data' filesep 'data_' country '_2022.mat']);
for t = 1:length(time_points)
    time_test = time_points(t);
    outdir = ['HELEN data1' filesep 'Epistatic networks' filesep country filesep char(time_test)];
    outfile = [outdir filesep 'edges.mat'];
    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    indSeq = (sampTimeSeqD <= time_test);
    M_t = M(indSeq,:);
    E = constructEpisNetwork(M_t,0.05);
    outStruct = struct(); 
    outStruct.E = E;
    parsave(outfile,outStruct);
end
%% 
% This section calculates and saves density-based p-values of VOCs and VOIs
% in temporal epistatic networks
% parpool(40);
nSubgrNode = 3000;
nMut = size(M,2);

outdirRes = ['HELEN data1' filesep 'p-values' filesep country];
if ~exist(outdirRes, 'dir')
    mkdir(outdirRes)
end
parfor t = 1:length(time_points)
    time_test = time_points(t)
    netData = load(['HELEN data1' filesep 'Epistatic networks' filesep country filesep char(time_test) filesep 'edges.mat']);
    E = netData.E;
    G = graph(E(:,1),E(:,2));
    if max(max(E)) < nMut
        G = addnode(G,nMut - max(max(E)));
    end
    for s = 1:size(voc,1)  
        outfile = [outdirRes filesep 'res_s' int2str(s) '_t' int2str(t) '.mat'];
        p_voc = densPValue(G,voc(s,:),nSubgrNode);
        stat_tp_v = struct(); 
        stat_tp_v.p_voc = p_voc;
        parsave(outfile,stat_tp_v);
    end
end
%% 
% This section find densest subgraphs in temporal epistatic networks, compare them
% with VOCs/VOIs, calculate the corresponding f-scores and saves the
% results

sens_densest = zeros(length(time_points),size(voc,1));
spec_densest = zeros(length(time_points),size(voc,1));
outdirRes = ['HELEN data1' filesep 'Densest subgraphs'];
if ~exist(outdirRes, 'dir')
    mkdir(outdirRes)
end
t = 0;
for time_test = time_points
    t = t+1
      
    netData = load(['HELEN data1' filesep 'Epistatic networks' filesep country filesep char(time_test) filesep 'edges.mat']);
    E = netData.E;        
    G = graph(E(:,1),E(:,2));
    if max(max(E)) < nMut
        G = addnode(G,nMut - max(max(E)));
    end
        
    subgr = getDensestSubgraph(G);
    subgr = find(subgr); 
    for i = 1:size(voc,1)
        sens_densest(t,i) = length(intersect(find(voc(i,:)),subgr))/sum(voc(i,:));
        spec_densest(t,i) = length(intersect(find(voc(i,:)),subgr))/length(subgr);
    end
end
f_densest = 2*sens_densest.*spec_densest./(sens_densest + spec_densest);
f_densest(isnan(f_densest)) = 0;
save([outdirRes filesep country '_densest_subgraph.mat'],'f_densest');

%% 
% This section infers viral variants with altered phenotypes in temporal epistatic networks using HELEN
% (Heralding Emerging Lineages in Epistatic Networks) framework

k_max = 27;
k_min = 7;
nSol = 100;
timeLimit = 20000;
stats = cell(1,length(time_points));
load(['HELEN data1' filesep 'Preprocessed input data' filesep 'data_' country '_2022.mat']);
outdirRes = ['HELEN data1' filesep 'Inferred variants' filesep country];
if ~exist(outdirRes, 'dir')
    mkdir(outdirRes)
end

for t = 1:length(time_points)
    time_test = time_points(t);
      
    netData = load(['HELEN data1' filesep 'Epistatic networks' filesep country filesep char(time_test) filesep 'edges.mat']);
    E = netData.E;        
    G = graph(E(:,1),E(:,2));
    if max(max(E)) < nMut
        G = addnode(G,nMut - max(max(E)));
    end
    outToken = int2str(t);
    [vocInfer, vocSupport] = HELEN_infer(G,k_min,k_max,nSol,timeLimit,outdirRes,outToken);
    
    allSens = zeros(length(vocInfer),size(voc,1));
    allSpec = zeros(length(vocInfer),size(voc,1));
    allF = zeros(length(vocInfer),size(voc,1));
    vocInferCummFreq = zeros(1,length(vocInfer));
    vocInferPrev = zeros(1,length(vocInfer));
    for j = 1:length(vocInfer)
        for i = 1:size(voc,1)
            allSens(j,i) = length(intersect(find(voc(i,:)),vocInfer{j}))/sum(voc(i,:));
            allSpec(j,i) = length(intersect(find(voc(i,:)),vocInfer{j}))/length(vocInfer{j});
            if allSens(j,i) + allSpec(j,i) ~= 0
                allF(j,i) = 2*allSens(j,i)*allSpec(j,i)/(allSens(j,i)+allSpec(j,i));
            else
                allF(j,i) = 0;
            end
        end
        ind1 = (sum(M(:,vocInfer{j}),2) == length(vocInfer{j}))';
        ind2 = (sampTimeSeqD <= time_points(t));
        vocInferCummFreq(j) = sum(ind1&ind2)/sum(ind2);
        if t >=2
            ind3 = (sampTimeSeqD > time_points(t-1));
            vocInferPrev(j) = sum(ind1&ind2&ind3)/sum(ind2&ind3);
        end
    end
    [valMax,indMax] = max(allF,[],2);
    stats{t} = struct();
    stats{t}.allF = allF;
    stats{t}.allSens = allSens;
    stats{t}.allSpec = allSpec;
    stats{t}.vocInfer = vocInfer;
    stats{t}.vocSupport = vocSupport;
    stats{t}.matchDist= vocSupport*valMax;
    stats{t}.vocInferCummFreq = vocInferCummFreq;
    stats{t}.vocInferPrev = vocInferPrev;
end
save(['HELEN data1' filesep 'Inferred variants' filesep country '_densestFixSizeExt.mat'],'stats');
