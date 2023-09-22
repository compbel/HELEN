function [varInfer, vocSupport] = HELEN_infer(G,k_min,k_max,nSol,timeLimit,outdirRes,outToken)

% generation of a pool of candidate dense communities of an epistatic network
nMut = numnodes(G);
commun2 = {};
[~,binsizes] = conncomp(G);
maxComp = max(binsizes);
k_min_curr = min([k_min maxComp]);
k = k_max;

while k >= k_min_curr
    if k > maxComp
        communCurr = [];
    else
        [solSets,~] = getPSubgraphExcl(G,k,commun2,nMut,nSol,timeLimit);
        communCurr = cell(1,size(solSets,1));
        for i = 1:size(solSets,1)
            communCurr{i} = find(solSets(i,:));
        end
    end
    if (~isempty(outdirRes)) && (~isempty(outToken))
        outfile = [outdirRes filesep 'commun_' outToken '_' int2str(k) '.mat'];
        outstruct = struct();
        outstruct.out = communCurr;
        parsave(outfile,outstruct);
    end
    commun2 = [commun2 communCurr];
    k = k-1
end

%assembly of viral variants from the pool of communities
[varInfer, vocSupport] = clustSubgraphs4(commun2,nMut,G);