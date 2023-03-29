function seqNew = fillGaps(seq,ref)
% nucs = 'ACTG';
nucs = 'ARNDCQEGHILKMFPSTWYV';
seqNew = seq;
for i = 1:size(seq,1)
    i
    for j = 1:size(seq,2)
        if ~contains(nucs,seq(i,j))
            seqNew(i,j) = ref(j);
        end
    end
end

