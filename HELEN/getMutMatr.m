function genotMatr = getMutMatr(hapAnalize,ref)

        n = size(hapAnalize,1);
        m = size(hapAnalize,2);
        genotMatr = zeros(n,m);
        for i = 1:n
            for j = 1:m
                if hapAnalize(i,j) ~= ref(j)
                    genotMatr(i,j) = 1;
                end
            end
        end