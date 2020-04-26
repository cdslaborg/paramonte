function CovMat = getCovMatFromCorMat(nd, StdVec, CorMat)

    CovMat = zeros(nd,nd);

    for j = 1 : nd
        CovMat(j,j) = StdVec(j)^2;
        for i = 1 : j-1
            CovMat(i,j) = CorMat(i,j) * StdVec(j) * StdVec(i);
            CovMat(j,i) = CovMat(i,j);
        end
    end

end