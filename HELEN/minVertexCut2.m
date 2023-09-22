function [optComps,optCut] = minVertexCut2(G,compThr,toDebug)

% Algorithm from "On computing the connectivities of graphs and digraphs"


    deg = degree(G);
    nodeIDs = 1:numnodes(G);
    while sum(deg == 1) > 0
        G = subgraph(G,deg > 1);
        nodeIDs = nodeIDs(deg > 1);
        deg = degree(G);
    end

    [bins,iC] = biconncomp(G,'OutputForm','cell');
    if ~isempty(iC)
        optCut = iC;
        optComps = bins(cellfun(@(c) length(c)>=compThr ,bins));
        for j = 1:length(optComps)
            optComps{j} = nodeIDs(optComps{j});
        end
        if toDebug
            figure
            h = plot(G,'layout','force');
            highlight(h,optCut,'NodeColor','r','MarkerSize',5)
            [];
        end
        return;
    end
    n = numnodes(G);
    deg = degree(G);
    A = adjacency(G);
    A2 = zeros(2*n,2*n);
    A2(1:n,(n+1):(2*n)) = eye(n,n);
    A2((n+1):(2*n),1:n) = A;
    labels = [1:n 1:n];
    
    s = find(deg==min(deg),1,'first');
    dist = distances(G,s);
    [dist,order] = sort(dist,'descend');
    L = find(dist>=3,1,'last');
    
%     AM_T = zeros(n,n);
%     AM_T(s,neighbors(G,s)) = 1;
%     leafs = false(1,n);
%     leafs(neighbors(G,s)) = 1;
%     treeNodes = false(1,n);
%     treeNodes([s neighbors(G,s)]) = 1;
%     while (sum(treeNodes) < n)
%     end
    
    
    
    optCut = [];
    optConnectivity = Inf;
    vertToCheck = ones(1,n);
    vertToCheck(s) = 0;
    for i = 1:L
        t = order(i);
        if (A(s,t) == 0) && (vertToCheck(t) == 1)
           disp([i sum(vertToCheck) optConnectivity])
           [connectivity_curr,cut_curr,GF] = minVertexCutPair(A2,s,t,n,labels);
           
            paths = zeros(connectivity_curr,n);
            for p = 1:connectivity_curr
                [VPe,dpe,EPe] = shortestpath(GF,s,t);
                P = unique(labels(VPe),'stable');
                paths(P(2:end-1)) = 1;
                GF = rmedge(GF,EPe);
            end
            ind = (prod(A*paths',2) > 0);
%             neighP = sum(A*paths' > 0,2);
%             ind = (neighP == connectivity_curr) | ((neighP == connectivity_curr-1)&A(:,s)) | ((neighP == connectivity_curr-1)&A(:,t));
            vertToCheck(ind) = 0;
            vertToCheck(t) = 0;
            ind1 = sum(A(:,ind),2) >= connectivity_curr;
            vertToCheck(ind1) = 0;
            
           
            if connectivity_curr < optConnectivity
                optConnectivity = connectivity_curr;
                optCut = cut_curr;
            end
            if optConnectivity == 2
                break;
            end
        end
    end
    
    if optConnectivity > 2
        Neigh = neighbors(G,s);
        for i = 1:length(Neigh)
            for j = (i+1):length(Neigh)
                s = Neigh(i);
                t = Neigh(j);
                if A(s,t) == 0
                    [i j]
                    [connectivity_curr,cut_curr] = minVertexCutPair(A2,s,t,n,labels);
                    if connectivity_curr < optConnectivity
                        optConnectivity = connectivity_curr;
                        optCut = cut_curr;
                    end
                    if optConnectivity == 2
                        break;
                    end
                end
            end
        end
    end
    
    if toDebug
        figure
        h = plot(G,'layout','force');
        highlight(h,optCut,'NodeColor','r','MarkerSize',5)
        [];
    end
    nodeIDsSub = setdiff(1:n,optCut,'stable');
    H = subgraph(G,nodeIDsSub);
    nodeIDs = nodeIDs(nodeIDsSub);
    bins = conncomp(H,'OutputForm','cell');
    
%     for v = optCut
%         for i = 1:length(bins)
%             if sum(A(v,nodeIDs(bins{i})),2) > 0
%                 bins{i} = union(bins{i},v);
%             end
%         end
%     end
        
    optComps = bins(cellfun(@(c) length(c)>=compThr ,bins));
    for j = 1:length(optComps)
        optComps{j} = nodeIDs(optComps{j});
    end
    

    function [connectivity,cut,GF] = minVertexCutPair(A2,s,t,n,labels)
        A2_curr = A2;
        A2_curr(s,:) = A2_curr(n+s,:);
        A2_curr(t,:) = A2_curr(n+t,:);
        ind = true(1,2*n);
        ind(n+s) = 0;
        ind(n+t) = 0;
        A2_curr = A2_curr(ind,ind);
        labels_curr = labels(ind);
        H = digraph(A2_curr);
        [mf,GF,cs,ct] = maxflow(H,s,t);
        [u,v] = find(A2_curr(cs,ct));
        u = cs(u);
        v = ct(v);
        connectivity = mf;
        cut = setdiff(unique([labels_curr(u) labels_curr(v)]),[s t]);