function subgraphs = randESU(G,k,sampProb,maxSampNode)
    global subgraphs;
    global countAll;
    global countNode;
    global A;
    global N;
    subgraphs = containers.Map('KeyType','int32','ValueType','any');
    countAll = int32(0);
    n = numnodes(G);
    N = cell(1,n);
    for i = 1:n
        N{i} = neighbors(G,i);
    end
    A = logical(adjacency(G));
    for i = 1:n
        countNode = 0;
        V_subgr = zeros(1,k);
        V_subgr(1) = i;
        currSize = 1;
%         V_ext = A(i,:);
%         V_ext(1:(i-1)) = 0;
        V_ext = N{i};
        V_ext = V_ext(V_ext > i);
        V_excl = union(N{i},1:i);
        extendSubgraph(V_subgr,currSize,V_ext,V_excl,k,sampProb,maxSampNode);
    end
    
function [] = extendSubgraph(V_subgr,currSize,V_ext,V_excl,k,sampProb,maxSampNode)
    global subgraphs;
    global countAll;
    global countNode;
    global A;
    global N;
    
    if countNode >= maxSampNode
        return;
    end
    
    if currSize == k
        countAll = countAll+1;
        countNode = countNode+1;
%         disp(countNode)
        subgraphs(countAll) = V_subgr;
        return;
    end
    while (~isempty(V_ext)) && (countNode < maxSampNode)
        if length(V_ext) > 1
            j = randsample(V_ext,1);
        else
            j = V_ext;
        end
        V_subgr_new = V_subgr;
        V_subgr_new(currSize+1) = j;
        V_ext = V_ext(V_ext~=j);
        V_ext_new = union(V_ext,setdiff(N{j},V_excl));
        V_excl_new = union(V_excl,N{j});
        p = rand;
        if p <= sampProb
            extendSubgraph(V_subgr_new,currSize+1,V_ext_new,V_excl_new,k,sampProb,maxSampNode);
        end
    end
    