function [optComps,nClust] = clustGraph(G,nClust,compThr,toDebug)

    f = figure('visible','off');
    h = plot(G,'Layout','force3');
    nClust = min([nClust size(h.XData,2)]);
    maxClust = min([2 size(h.XData,2)]);
    [clusters,~,D] = spectralcluster([h.XData' h.YData' h.ZData'],nClust);    
    gaps = abs(D(2:end)./D(1:end-1));
    nClust = find(gaps==max(gaps));
    [clusters,~,~] = spectralcluster([h.XData' h.YData' h.ZData'],nClust);

%     AM = adjacency(G);
%     [clusters,~,D] = spectralcluster(AM,10,'Distance','precomputed');
%     gaps = abs(D(2:end)./D(1:end-1));
%     nClust = find(gaps==max(gaps));
%     [clusters,~,D] = spectralcluster(AM,nClust,'Distance','precomputed');

    if nClust == 1
        myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton','replicate',10));
        eva = evalclusters([h.XData' h.YData' h.ZData'],myfunc,'gap','KList',1:maxClust);
        if (eva.OptimalK > 1) && (sum(isinf(eva.CriterionValues))==0)
            nClust = eva.OptimalK;
            clusters = eva.OptimalY;
        end
    end
    
%     
%     if nClust == 1
%         D = distances(G);
%         Y = cmdscale(D,4);
%         [clusters,~,D] = spectralcluster(Y,min([100 size(Y,1)]));    
%         gaps = abs(D(2:end)./D(1:end-1));
%         nClust = find(gaps==max(gaps));
%         [clusters,~,~] = spectralcluster(Y,nClust);
%     end


    if toDebug
        c = rand(nClust,3);
        colors = zeros(length(h.XData),3);
        for ii = 1:length(h.XData)
            colors(ii,:) = c(clusters(ii),:);
        end
        figure 
        scatter3(h.XData',h.YData',h.ZData',30,colors,'filled');
        [];
    end
    
    optComps = cell(1,nClust);
    for i = 1:nClust
        optComps{i} = find(clusters == i);
    end
    optComps = optComps(cellfun(@(c) length(c)>=compThr ,optComps));
    [];

 
        

  
