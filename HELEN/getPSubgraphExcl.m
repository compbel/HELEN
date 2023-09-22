function [solSets,solNEdges] = getPSubgraphExcl(G,k,exclSet,nMut,nSol,timeLimit)
n = numnodes(G);
E = G.Edges.EndNodes;
m = size(E,1);

id = randi(1000000);
filename = ['ILP_' int2str(k) '_' int2str(id) '.lp'];
fileID = fopen(filename,'w');
fprintf(fileID,'%s\n','maximize');
for e = 1:m
    Yij = ['Y' int2str(E(e,1)) ',' int2str(E(e,2))];
    fprintf(fileID,'%s',['+ ' Yij ' ']);
end
fprintf(fileID,'\n%s\n',' st');
for i = 1:n
    Xi = ['X' int2str(i)];
    fprintf(fileID,'%s',['+ ' Xi ' ']);     
end
fprintf(fileID,'%s\n',['= ' int2str(k)]);

for e = 1:m
    Yij = ['Y' int2str(E(e,1)) ',' int2str(E(e,2))];
    Xi = ['X' int2str(E(e,1))];
    Xj = ['X' int2str(E(e,2))];
    fprintf(fileID,'%s\n',[Yij ' - ' Xi ' <= 0']);
    fprintf(fileID,'%s\n',[Yij ' - ' Xj ' <= 0']);
end

for i = 1:length(exclSet)
    set = setdiff(1:nMut,exclSet{i});
    for j = set
        Xj = ['X' int2str(j)];
        fprintf(fileID,'%s',['+ ' Xj ' ']);  
    end
    fprintf(fileID,'%s\n',['>= 1']);
end
fprintf(fileID,'%s\n','binaries');
for i = 1:n
    Xi = ['X' int2str(i)];
    fprintf(fileID,'%s\n',Xi);
end
for e = 1:m
    Yij = ['Y' int2str(E(e,1)) ',' int2str(E(e,2))];
    fprintf(fileID,'%s\n',Yij);
end

fprintf(fileID,'%s\n','end');
fclose(fileID);

clear model;
model = gurobi_read(filename);

params.outputflag = 1;
params.timelimit = timeLimit;
% params.NodefileStart = 30;
% params.MIPFocus = 2;
params.threads = 64;
params.PoolSearchMode = 2;
params.PoolSolutions = nSol;
% params.PoolGap = 0.4;
result = gurobi(model, params);
if ~strcmp(result.status,'INFEASIBLE')
    solSets = false(n,size(result.pool,2));
    solNEdges = zeros(1,size(result.pool,2)); 
    for i = 1:size(result.pool,2)
        solSets(:,i) = round(result.pool(i).xn((m+1):end));
        solNEdges(i) = result.pool(i).objval;
    end
    [solSets,ia,ic] = unique(solSets','rows');
    solNEdges = solNEdges(ia);
else
    solSets = [];
    solNEdges = [];
end
delete(filename);
