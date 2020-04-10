% This subroutine is the same as combineCovMean, with the IMPORTANT difference that
% only the upper triangles and diagonals of the input covariance matrices need to be given by the user: CovMatUpper, CovMatUpperB
function [MeanVec, CovMatUpper] = combineMeanCovUpper(nd, npA, MeanVecA, CovMatUpperA, npB, MeanVecB, CovMatUpperB)

    npAreal     = npA;% xxx perhaps remove these conversions, may not have significant performance benefit
    npBreal     = npB;% xxx perhaps remove these conversions, may not have significant performance benefit
    npABinverse = 1 / (npA + npB);

    % First find MeanVec:
    MeanVec     = npABinverse * (npAreal*MeanVecA + npBreal*MeanVecB);
    
    % Now find new Covariance matrix:
    for j = 1 : nd
        for i = 1 : j
            CovMatUpper(i,j)  = ( npAreal * ( CovMatUpperA(i,j) + MeanVecA(i) * MeanVecA(j) )   ...
                                + npBreal * ( CovMatUpperB(i,j) + MeanVecB(i) * MeanVecB(j) )   ...
                                ) * npABinverse - MeanVec(i) * MeanVec(j)                       ...
                                ;
        end
    end

end