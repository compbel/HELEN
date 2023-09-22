function p_voc = densPValue(G,voc,nSubgrNode)

nMut = numnodes(G);
nonisol = find(degree(G) > 0);
G1 = subgraph(G,nonisol);
labels = arrayfun(@num2str, 1:nMut, 'UniformOutput', 0);
mut_G1 = cellfun(@str2num,labels(nonisol));

[~,idx]=ismember(find(voc),nonisol);
H = subgraph(G1,idx(idx>0));
labels_H = labels(nonisol(idx(idx > 0)));
        
if ~isempty(labels_H)
    [bins,binsizes] = conncomp(H);
    ind = find(binsizes == max(binsizes));
    H = subgraph(H,bins == ind(1));
    mut_H = cellfun(@str2num,labels_H(bins == ind(1)));

    densStr = numedges(H)/numnodes(H);
    nDens = 0;
    subgraphs = randESU(G1,numnodes(H),1,nSubgrNode);
    AM1 = adjacency(G1);
    for j = 1:size(subgraphs,1)
        AM_samp = AM1(subgraphs(j),subgraphs(j));
        if sum(sum(AM_samp))/(2*size(AM_samp,1)) >= densStr
            if ~isequal(sort(mut_G1(subgraphs(j))),mut_H)
                nDens = nDens + 1;
            end
        end
    end
    p_voc = nDens/size(subgraphs,1);
else
    p_voc = 1;
end