% This subroutine is exact same thing as getSamCovUpperMeanTrans,
% with the only difference that here points can have different weights. This subroutine is specifically optimized for use in ParaMCMC.
function [CovMatUpper, Mean] = getWeiSamCovUppMeanTrans(np, sumWeight, nd, Point, Weight)

    % np,nd,sumWeight        : np is the number of observations, nd is the number of parameters for each observation
    % Weight(np)             : weights of the points
    % Point(nd,np)           : Point is the matrix of the data, CovMatUpper contains the elements of the sample covariance matrix
    % CovMatUpper(nd,nd)     : Covariance matrix of the input data
    % Mean(nd)               : Mean vector

    CovMatUpper = zeros(nd,nd);
    NormedData  = zeros(nd,np);
    Mean        = zeros(nd, 1);

    for i = 1 : np
        for j = 1 : nd
            Mean(j) = Mean(j) + Weight(i) * Point(j,i);
        end
    end
    Mean = Mean / sumWeight;

    for i = 1 : np
        NormedData(1:nd,i) = Point(1:nd,i) - Mean;
    end

    sumWeightMinusOneInvReal = 1 / (sumWeight-1);
    for j = 1 : nd
        for i = 1 : j
            CovMatUpper(i,j) = 0;
            for ip = 1 : np
                CovMatUpper(i,j) = CovMatUpper(i,j) + Weight(ip) * NormedData(i,ip) * NormedData(j,ip);
            end
            CovMatUpper(i,j) = CovMatUpper(i,j) * sumWeightMinusOneInvReal;
        end
    end

    Mean = Mean';

end
