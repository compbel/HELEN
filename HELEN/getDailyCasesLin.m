function dailyCasesLin = getDailyCasesLin(timeObsCases,sampTimeSeq)

dailyCasesLin = zeros(length(timeObsCases),1);
for t = 1:length(timeObsCases)
    dailyCasesLin(t) = sum(sampTimeSeq == t);
end