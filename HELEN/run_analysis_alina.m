clear;

% path to GUROBI LP solver
addpath 'C:\gurobi811\win64\matlab'
addpath '/opt/gurobi810/linux64/matlab'


continents = {'Europe','Europe','North_America','North_America','Asia','Asia','Africa','Africa','Oceania','South_America','South_America','Europe','Europe','Europe','Europe','Europe'};
countries = {'United_Kingdom','Germany','USA','Canada','Japan','India','South_Africa','Kenya','Australia','Brazil','Peru','France','Sweden','Spain','Denmark','Italy'};

time_test_begin = datetime('05-01-2020','InputFormat','MM-dd-yyyy');
time_test_end = datetime('12-31-2021','InputFormat','MM-dd-yyyy');
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
outfolder = 'HELEN_data1';
parpool(45);

voc_names = {'Alpha','Beta','Gamma', 'Delta','Omicron','Lambda','Mu','Theta','Eta','Kappa','Zeta','Epsilon','Iota'};
firstSamp = {'01-Sep-2020','01-May-2020','01-Nov-2020','01-Oct-2020','01-Nov-2021','01-Aug-2020','01-Jan-2021','01-Jan-2021','01-Dec-2020','01-Oct-2020','17-Aug-2021','09-Nov-2021','20-Sep-2021'};
firstSamp = datetime(firstSamp,'InputFormat','dd-MMM-yyyy');
% 
for cnt = 1:length(countries)
    % Set the country to be analyzed
    continent = continents{cnt};
    country = countries{cnt}

    % This section loads the data, prepares several data structures used by
    % HELEN scripts and saves them for further use

    % Load sequencing data and metadata downloaded from GISAID
    outdir = [outfolder filesep 'Preprocessed_input_data'];
    if ~exist([outfolder filesep 'Preprocessed_input_data'], 'dir')
        mkdir([outfolder filesep 'Preprocessed_input_data'])
    end
    if ~exist([outdir filesep 'data_' country '_2022.mat'],'file')
        seqF = fastaread(['Genomic_data_upd' filesep continent filesep 'sequences_' country '_protein.fasta']);
        seq = char(seqF.Sequence);
        seq = seq(:,1:end-1);

        % Parse sampling times of sequences from GISAID data
        sampTimeSeqD = NaT(1,length(seqF));
        sampTimeSeq = zeros(1,length(seqF));
        epiID = zeros(1,length(seqF));
        strain = zeros(1,length(seqF));
        badDate = false(1,length(seqF));
        formatDate = 'yyyy-MM-dd';
        parfor i = 1:length(seqF)
            i
            C = strsplit(seqF(i).Header,'|');
            dateSeq = C{3};
            epiID(i) = str2double(C{2}(9:end));
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
        epiID = epiID(ind);

        ind = (sampTimeSeqD <= time_points(end));
        seqF = seqF(ind);
        seq = seq(ind,:);
        sampTimeSeq = sampTimeSeq(ind);
        sampTimeSeqD = sampTimeSeqD(ind);
        epiID = epiID(ind);

        % construct mutation matrix
        refF = fastaread(['Genomic_data_upd' filesep 'reference_s_aa.fas']);
        ref = char(refF.Sequence);
        seq = fillGaps(seq,ref);
        M = getMutMatr(seq,ref);
        M = sparse(M);

        % clear memory and save 
        if exist([outfolder filesep 'Preprocessed_input_data' filesep 'sequences_' country '_2022.fas'], 'file')
            delete([outfolder filesep 'Preprocessed_input_data' filesep 'sequences_' country '_2022.fas']);
        end
        fastawrite([outfolder filesep 'Preprocessed_input_data' filesep 'sequences_' country '_2022.fas'],seqF);
        clear seq seqF aux badDate C dataEpid dateSeq ind sampCounts 
        if ~exist([outfolder filesep 'Preprocessed_input_data'], 'dir')
            mkdir([outfolder filesep 'Preprocessed_input_data'])
        end
        outStruct = struct();
        outStruct.M = M;
        outStruct.sampTimeSeqD = sampTimeSeqD;
        outStruct.epiID = epiID;
        save([outfolder filesep 'Preprocessed_input_data' filesep 'data_' country '_2022.mat'],'-v7.3');
        parsaveXL([outfolder filesep 'Preprocessed_input_data' filesep 'data_' country '_2022.mat'],outStruct);
    end
    metadataCont = readtable(['Genomic_data_upd' filesep 'Metadata' filesep 'metadata_2022_05_01_updated_' continent '_ui_extracted.tsv'],"FileType","text",'Delimiter', '\t');
    ind = strcmp(metadataCont.Country,country);
    metadata = metadataCont(ind,:);

    pat = digitsPattern;
    epiIDAll = extract(metadata.AccessionID,pat);
    epiIDAll = str2double(epiIDAll);
    underInvestAll = metadata.UnderInvestigation;

    strainAll = zeros(1,length(epiIDAll));
    for i = 1:length(voc_names)
        i
        ind = strcmp(metadata.specifiedVariant,voc_names{i});
        ind = (metadata.categoricalVariants == i);
        strainAll(ind) = i;
    end
    save(['Preprocessed_input_data' filesep 'metadata_' country '.mat'],'epiIDAll','strainAll','underInvestAll');

    underInvest = false(1,size(M,1));
    suspicious = false(1,size(M,1));
    strainGIS = zeros(1,size(M,1));
    found = true(1,size(M,1));
    parfor i = 1:size(M,1)
        i
        j = find(epiIDAll == epiID(i));
        if isempty(j)
            continue;
        end
        found(i) = 1;
        strainGIS(i) = strainAll(j);
        underInvest(i) = underInvestAll(j);
        if (strainGIS(i)>0) && (strainGIS(i)<=5) && (sampTimeSeqD(i) < firstSamp(strainGIS(i)))
            suspicious(i) = 1;
        end
    end
    save([outfolder filesep 'Preprocessed_input_data' filesep 'seq_suspect_' country '_2022.mat'],'underInvest','suspicious');
end


%% 

% clear;
% 
% %path to GUROBI LP solver
% % addpath 'C:\gurobi811\win64\matlab'
% addpath '/opt/gurobi810/linux64/matlab'
% 
% 
% continents = {'Europe','Europe','North_America','North_America','Asia','Asia','Africa','Africa','Oceania','South_America','South_America','Europe','Europe','Europe','Europe','Europe'};
% countries = {'United_Kingdom','Germany','USA','Canada','Japan','India','South_Africa','Kenya','Australia','Brazil','Peru','France','Sweden','Spain','Denmark','Italy'};
% 
% time_test_begin = datetime('05-01-2020','InputFormat','MM-dd-yyyy');
% time_test_end = datetime('12-31-2021','InputFormat','MM-dd-yyyy');
% time_points = time_test_begin:caldays(15):time_test_end;
% 
% % voc/voi haplotypes
% mut_alpha = [501 570 614 681 716 982 1118];
% mut_beta = [80 215 417 484 501 614 701];
% mut_gamma = [18 20 26 138 190 417 484 501 614 655 1027 1176];
% mut_delta = [19 142 158 452 478 614 681 950];
% mut_omicron = [67 95 145 212 339 371 373 375 417 440 446 477 478 484 493 496 498 501 505 547 614 655 679 681 764 796 856 954 969 981];
% mut_lambda = [75 76 452 490 614 859];
% mut_mu = [95 145 346 484 501 614 681 950];
% mut_theta = [265 484 501 614 681 1092 1101 1176];
% mut_eta = [52 67 484 614 677 888];
% mut_kappa = [142 154 452 484 614 681 1071];
% mut_zeta = [484 614 1176];
% mut_epsilon = [13 152 452 614];
% mut_iota = [5 95 614];
% 
% nMut = 1273;
% voc = false(10,nMut);
% voc(1,mut_alpha) = 1;
% voc(2,mut_beta) = 1;
% voc(3,mut_gamma) = 1;
% voc(4,mut_delta) = 1;
% voc(5,mut_omicron) = 1;
% voc(6,mut_lambda) = 1;
% voc(7,mut_mu) = 1;
% voc(8,mut_theta) = 1;
% voc(9,mut_eta) = 1;
% voc(10,mut_kappa) = 1;
% 
% % days of first reported COVID-19 cases
% day1 = {'30-Jan-2020','26-Jan-2020','21-Jan-2020','16-Feb-2020','21-Jan-2020','29-Jan-2020','04-Mar-2020','12-Mar-2020',...
%     '23-Feb-2020','25-Feb-2020','05-Mar-2020','23-Jan-2020','31-Jan-2020','31-Jan-2020','26-Feb-2020','30-Jan-2020'};
% day1 = datetime(day1,'InputFormat','dd-MMM-yyyy');
% 
% outfolder = 'HELEN_data1';
% parpool(45);
% % 
% 
% for cnt = 1:length(countries)
%     country = countries{cnt};
%     % This section constructs epistatic networks and saves them for further use
%     load([outfolder filesep 'Preprocessed_input_data' filesep 'data_' country '_2022.mat']);
% %     load([outfolder filesep 'Preprocessed_input_data' filesep 'seq_suspect_' country '_2022.mat']);
% %     ind = (~suspicious)&(~underInvest);
% % %     ind = (~underInvest);
% %     M = M(ind,:);
% %     sampTimeSeqD = sampTimeSeqD(ind);
%     
%     for t = 1:length(time_points)
%         time_test = time_points(t);
%         outdir = [outfolder filesep 'Epistatic_networks' filesep country filesep char(time_test)];
%         outfile = [outdir filesep 'edges.mat'];
%         if ~exist(outfile,'file')
%             if ~exist(outdir, 'dir')
%                 mkdir(outdir)
%             end
%             indSeq = (sampTimeSeqD <= time_test);
%             M_t = M(indSeq,:);
%             E = constructEpisNetwork(M_t,0.05);
%             outStruct = struct(); 
%             outStruct.E = E;
%             parsave(outfile,outStruct);
%         end
%     end
% end
%% 


% [X,Y,Z] = meshgrid(1:length(countries),1:length(time_points),1:size(voc,1));
% triples = [X(:) Y(:) Z(:)];
% 
% parfor tr = 1:length(triples)
%     cnt = triples(tr,1);
%     t = triples(tr,2);
%     s = triples(tr,3);
%     country = countries{cnt};
%     [cnt s t]
%     % This section calculates and saves density-based p-values of VOCs and VOIs
%     % in temporal epistatic networks
%     nSubgrNode = 3000;
% 
%     outdirRes = [outfolder filesep 'p-values' filesep country];
%     if ~exist(outdirRes, 'dir')
%         mkdir(outdirRes)
%     end    
%     time_test = time_points(t);
%     netData = load([outfolder filesep 'Epistatic_networks' filesep country filesep char(time_test) filesep 'edges.mat']);
%     E = netData.E;
%     G = graph(E(:,1),E(:,2));
%     if max(max(E)) < nMut
%         G = addnode(G,nMut - max(max(E)));
%     end
%     outfile = [outdirRes filesep 'res_s' int2str(s) '_t' int2str(t) '.mat'];
%     if ~exist(outfile,'file')
%         disp(['processing ' num2str([cnt s t])])
%         p_voc = densPValue(G,voc(s,:),nSubgrNode);
%         stat_tp_v = struct(); 
%         stat_tp_v.p_voc = p_voc;
%         parsave(outfile,stat_tp_v);
%     end
% end
%% 

% 
% 
% for cnt = 1:length(countries)
%     % This section find densest subgraphs in temporal epistatic networks, compare them
%     % with VOCs/VOIs, calculate the corresponding f-scores and saves the
%     % results
% 
%     country = countries{cnt};
%     sens_densest = zeros(length(time_points),size(voc,1));
%     spec_densest = zeros(length(time_points),size(voc,1));
%     outdirRes = [outfolder filesep 'Densest_subgraphs'];
%     if ~exist(outdirRes, 'dir')
%         mkdir(outdirRes)
%     end
%     t = 0;
%     for time_test = time_points
%         t = t+1
% 
%         netData = load([outfolder filesep 'Epistatic_networks' filesep country filesep char(time_test) filesep 'edges.mat']);
%         E = netData.E;        
%         G = graph(E(:,1),E(:,2));
%         if max(max(E)) < nMut
%             G = addnode(G,nMut - max(max(E)));
%         end
% 
%         subgr = getDensestSubgraph(G);
%         subgr = find(subgr); 
%         for i = 1:size(voc,1)
%             sens_densest(t,i) = length(intersect(find(voc(i,:)),subgr))/sum(voc(i,:));
%             spec_densest(t,i) = length(intersect(find(voc(i,:)),subgr))/length(subgr);
%         end
%     end
%     f_densest = 2*sens_densest.*spec_densest./(sens_densest + spec_densest);
%     f_densest(isnan(f_densest)) = 0;
%     save([outdirRes filesep country '_densest_subgraph.mat'],'f_densest');
% end
%% 

clear;
%path to GUROBI LP solver
% addpath 'C:\gurobi811\win64\matlab'
addpath '/opt/gurobi810/linux64/matlab'


continents = {'Europe','Europe','North_America','North_America','Asia','Asia','Africa','Africa','Oceania','South_America','South_America','Europe','Europe','Europe','Europe','Europe'};
countries = {'United_Kingdom','Germany','USA','Canada','Japan','India','South_Africa','Kenya','Australia','Brazil','Peru','France','Sweden','Spain','Denmark','Italy'};

time_test_begin = datetime('05-01-2020','InputFormat','MM-dd-yyyy');
time_test_end = datetime('12-31-2021','InputFormat','MM-dd-yyyy');
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
% parpool(45);
outfolder = 'HELEN_data3';
%% 


[X,Y] = meshgrid(1:length(countries),1:length(time_points));
pairs = [X(:) Y(:)];

parfor pr = 1:length(pairs)
    cnt = pairs(pr,1);
    country = countries{cnt};
    t = pairs(pr,2);

    % This section infers viral variants with altered phenotypes in temporal epistatic networks using HELEN
    % (Heralding Emerging Lineages in Epistatic Networks) framework

    k_max = 27;
    k_min = 7;
    nSol = 100;
    timeLimit = 20000;
    x = load([outfolder filesep 'Preprocessed_input_data' filesep 'data_' country '_2022.mat']);
    sampTimeSeqD = x.sampTimeSeqD;
    M = x.M;
    outdirRes = [outfolder filesep 'Inferred_variants' filesep country];
    if ~exist(outdirRes, 'dir')
        mkdir(outdirRes)
    end

        time_test = time_points(t);
        netData = load([outfolder filesep 'Epistatic_networks' filesep country filesep char(time_test) filesep 'edges.mat']);
        E = netData.E;        
        G = graph(E(:,1),E(:,2));
        if max(max(E)) < nMut
            G = addnode(G,nMut - max(max(E)));
        end
        outToken = int2str(t);
        
        outfile = [outdirRes filesep 'commun_' outToken '.mat'];
        if ~exist(outfile, 'file')

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
            stat = struct();
            stat.allF = allF;
            stat.allSens = allSens;
            stat.allSpec = allSpec;
            stat.vocInfer = vocInfer;
            stat.vocSupport = vocSupport;
            stat.matchDist= vocSupport*valMax;
            stat.vocInferCummFreq = vocInferCummFreq;
            stat.vocInferPrev = vocInferPrev;
            parsave(outfile,stat);
        end
end

%% 

% 
% 
% for cnt = 1:length(countries)
%     country = countries{cnt}
%     outdirRes = [outfolder filesep 'Inferred_variants' filesep country];
%     stats = cell(1,length(time_points));
%     for t = 1:length(time_points)
%         outToken = int2str(t);
%         load([outdirRes filesep 'commun_' outToken '.mat']);
%         stats{t} = struct();
%         stats{t}.allF = allF;
%         stats{t}.allSens = allSens;
%         stats{t}.allSpec = allSpec;
%         stats{t}.vocInfer = vocInfer;
%         stats{t}.vocSupport = vocSupport;
%         stats{t}.matchDist= matchDist;
%         stats{t}.vocInferCummFreq = vocInferCummFreq;
%         stats{t}.vocInferPrev = vocInferPrev;
%     end
%     save([outfolder filesep 'Inferred_variants' filesep country '_densestFixSizeExt.mat'],'stats');
% end
