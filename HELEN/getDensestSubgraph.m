function subgr = getDensestSubgraph(G)
n = numnodes(G);
E = G.Edges.EndNodes;
m = size(E,1);

filename = ['ILP_densest' num2str(randi(1000000)) '.lp'];
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
fprintf(fileID,'%s\n',['<= 1']);

for e = 1:m
    Yij = ['Y' int2str(E(e,1)) ',' int2str(E(e,2))];
    Xi = ['X' int2str(E(e,1))];
    Xj = ['X' int2str(E(e,2))];
    fprintf(fileID,'%s\n',[Yij ' - ' Xi ' <= 0']);
    fprintf(fileID,'%s\n',[Yij ' - ' Xj ' <= 0']);
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
model = gurobi_relax(model);

params.outputflag = 1;
params.timelimit = 20000;
% params.NodefileStart = 30;
% params.MIPFocus = 2;
params.threads = 18;
% params.PoolSearchMode = 2;
% params.PoolSolutions = 1000;
% params.PoolGap = 0.4;
result = gurobi(model, params);
subgr = result.x(end-n+1:end) > 0;

delete(filename); 
[];
