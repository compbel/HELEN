function [vocInfer,vocFreq] = clustSubgraphs4(subgraphs,nMut,G)
compSizeThr = 1;
compSizeThrCut = 5;
lenThr = 5;
nClust = 10;
toDebug = false;

    if isa(subgraphs,'containers.Map')
        subgraphs = values(subgraphs);
    end
    subgrsV = zeros(length(subgraphs),nMut);
    for i = 1:length(subgraphs)
        G_sub = subgraph(G,subgraphs{i});
        [comps,compsizes] = conncomp(G_sub);
        lgComp = find(compsizes == max(compsizes),1,'first');
        subgraphs{i} = subgraphs{i}(comps == lgComp);
        subgrsV(i,subgraphs{i}) = 1;
    end   
    [subgrsV,ia,~] = unique(subgrsV,'rows');
    subgraphs = subgraphs(ia);
    
    subgrsVNew = zeros(max(sum(subgrsV,2))*length(subgraphs),nMut);
    subgraphsNew = cell(1,max(sum(subgrsV,2))*length(subgraphs));
    nSub = 0;
    for i = 1:length(subgraphs)
        G_sub = subgraph(G,subgraphs{i});
        [bins,~] = biconncomp(G_sub,'OutputForm','cell');
        bigBins = bins(cellfun(@(c) length(c)>=compSizeThrCut ,bins));
        maxBinSize = max(cellfun(@(c) length(c) ,bins));
        if (~isempty(bigBins)) && ((maxBinSize~=length(subgraphs{i})-1) || length(subgraphs)>=2*compSizeThrCut)
            for j = 1:length(bigBins)
                nSub = nSub + 1;
                subgraphsNew{nSub} = subgraphs{i}(bigBins{j});
                subgrsVNew(nSub,subgraphsNew{nSub}) = 1;
            end
        else
            nSub = nSub + 1;
            subgraphsNew{nSub} = subgraphs{i};
            subgrsVNew(nSub,subgraphsNew{nSub}) = 1;
        end
    end
    subgrsV = subgrsVNew(1:nSub,:);
    subgraphs = subgraphsNew(1:nSub);
    [subgrsV,ia,~] = unique(subgrsV,'rows');
    subgraphs = subgraphs(ia);
    
    DM = subgrsV*subgrsV';
    
    S = repmat(sum(subgrsV,2),1,length(DM));
    
    AM = (DM == (min(S,S')-1));
    H = graph(AM);

    
    [comps,compsizes] = conncomp(H,'OutputForm','cell');
    comps = comps(compsizes>=compSizeThr);
    compsizes = compsizes(compsizes>=compSizeThr);
    
    if toDebug
        [];
    end
    
    
    for i = 1:length(compsizes)
        if compsizes(i) >= 2*compSizeThrCut
            comp = comps{i};
            [optComps,~] = minVertexCut2(subgraph(H,comp),compSizeThrCut,toDebug);
            ind = [i (length(comps)+1):(length(comps)+length(optComps)-1)];
            for j = 1:length(optComps)
                comps{ind(j)} = comp(optComps{j});
                compsizes(ind(j)) = length(optComps{j});
            end
        end
    end
    
    
    toClust = true(1,length(compsizes));
    for iter = 1:3
        for i = 1:length(compsizes)
            if (compsizes(i) >= 2*compSizeThrCut) && (toClust(i))
                comp = comps{i};
                [optComps,nClustComp] = clustGraph(subgraph(H,comp),nClust,compSizeThrCut,toDebug);
                ind = [i (length(comps)+1):(length(comps)+length(optComps)-1)];
                for j = 1:length(optComps)
                    comps{ind(j)} = comp(optComps{j});
                    compsizes(ind(j)) = length(optComps{j});
                    if nClustComp > 1
                        toClust(ind(j)) = 1;
                    else
                        toClust(ind(j)) = 0;
                    end
                end
            else
                toClust(i) = 0;
            end
        end
    end
            
    compSizeThr = min([compSizeThr max(compsizes)]);
    vocInfer = cell(1,sum(compsizes >= compSizeThr));
    vocFreq = zeros(1,sum(compsizes >= compSizeThr));
    for i = 1:length(compsizes)        
        comp = comps{i};
        vocFreq(i) = length(comp);
        muts = find(sum(subgrsV(comp,:),1)>=min([2 max(sum(subgrsV(comp,:),1))]));
        if length(muts) == 1
            muts = find(sum(subgrsV(comp,:),1)>0);
        end

        
        if length(comp) == 1
            vocInfer{i} = muts;
            continue;
        end
        
        G_comp = subgraph(G,muts);
        deg = degree(G_comp);
        nClust = Inf;
        for i1 = 1:5
            eva = evalclusters(deg,'linkage','gap','KList',1:3);
            if sum(isinf(eva.CriterionValues)) > 0
                continue;
            end
            if eva.OptimalK < nClust
                nClust = eva.OptimalK;
            end
        end
        if isinf(nClust)
            nClust = 1;
        end
        
        while nClust >= 1
            [clust,centr] = kmeans(deg,nClust,'Replicates',10);
            lgClust = find(centr == max(centr));
            if sum(clust == lgClust) > 1
                break;
            else
                nClust = nClust - 1;
            end
        end
        ind = (clust == lgClust);
        G_comp = subgraph(G_comp,ind);
        muts = muts(ind);
                
        subgr = getDensestSubgraph(G_comp);
        vocInfer{i} = muts(subgr);             
    end
    ind = cellfun(@isempty,vocInfer);
    vocInfer(ind) = [];
    vocFreq(ind) = [];
    
    [~,ia,ic] = unique(cellfun(@num2str,vocInfer,'uni',0),'stable');
    vocInfer = vocInfer(ia);
    vocFreq1 = zeros(1,length(ia));
    for i = 1:length(ia)
        vocFreq1(i) = sum(vocFreq(ic==i));
    end
    
    ind = (cellfun('length',vocInfer)>=lenThr);
    vocInfer = vocInfer(ind);
    vocFreq1 = vocFreq1(ind);
    vocFreq = vocFreq1; 
    vocFreq = vocFreq/sum(vocFreq);    
    