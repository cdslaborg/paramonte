% This subroutine is exact same thing as getSamCovMeanTrans, with the only differences that only the upper triangle of the covariance matrix is returned.
% Also, optional arguments are not available. This subroutine is specifically optimized for use in ParaMCMC.

function [CovMatUpper, Mean] = getSamCovUpperMeanTrans(np,nd,Point)

    % np,nd                  : np is the number of observations, nd is the number of parameters for each observation
    % Point(nd,np)           : Point is the matrix of the data, CovMatUpper contains the elements of the sample covariance matrix
    % CovMatUpper(nd,nd)     : Covariance matrix of the input data
    % Mean(nd)               : Mean vector

    CovMatUpper = zeros(nd,nd);
    NormedData  = zeros(nd,np);
    Mean        = zeros(nd, 1);

    for i = 1 : np
        for j = 1 : nd
            Mean(j) = Mean(j) + Point(j,i);
        end
    end
    Mean = Mean / np;

    for i = 1 : np
        NormedData(1:nd,i) = Point(1:nd,i) - Mean;
    end

    npMinusOneInvReal = 1 / (np-1);
    for j = 1 : nd
        for i = 1 : j
            CovMatUpper(i,j) = dot(NormedData(i,1:np), NormedData(j,1:np)) * npMinusOneInvReal;
        end
    end

    Mean = Mean';

end