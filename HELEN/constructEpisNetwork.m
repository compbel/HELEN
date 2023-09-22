function E = constructEpisNetwork(M,rho)

[X,Y] = meshgrid(1:size(M,2),1:size(M,2));
pairs = [X(:) Y(:)];
episPairs = false(1,length(pairs));

parfor pair = 1:length(pairs)
    if pairs(pair,1) >= pairs(pair,2)
        continue;
    end
    m1 = pairs(pair,1);
    m2 = pairs(pair,2);
    [m1 m2]
    lineage00_pair = logical((1-M(:,m1)).*(1-M(:,m2)));
    lineage10_pair = (M(:,m1) - M(:,m2)) > 0;
    lineage01_pair = (M(:,m2) - M(:,m1)) > 0;
    lineage11_pair = logical(M(:,m1).*M(:,m2));

    O_00 = sum(lineage00_pair); 
    O_10 = sum(lineage10_pair); 
    O_01 = sum(lineage01_pair); 
    O_11 = sum(lineage11_pair); 
    n = size(M,1);
    L = size(M,2);
    p = (O_10*O_01)/(O_00*n);
    prob = 1 - binocdf(O_11-1,n,p);
    thr = rho/nchoosek(L,2);
    if (prob <= thr) && (O_11 >= 1) && (O_10 >= 1) && (O_01 >= 1)
        episPairs(pair) = 1;       
    end
end
E = pairs(episPairs,:);

