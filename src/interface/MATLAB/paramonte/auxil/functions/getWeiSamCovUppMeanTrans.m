%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ParaMonte: plain powerful parallel Monte Carlo library.
%
%   Copyright (C) 2012-present, The Computational Data Science Lab
%
%   This file is part of the ParaMonte library.
%
%   ParaMonte is free software: you can redistribute it and/or modify it 
%   under the terms of the GNU Lesser General Public License as published 
%   by the Free Software Foundation, version 3 of the License.
%
%   ParaMonte is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with the ParaMonte library. If not, see, 
%
%       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
%
%   ACKNOWLEDGMENT
%
%   As per the ParaMonte library license agreement terms, 
%   if you use any parts of this library for any purposes, 
%   we ask you to acknowledge the ParaMonte library's usage
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This subroutine is exact same thing as getSamCovUpperMeanTrans,
% with the only difference that here points can have different weights. This subroutine is specifically optimized for use in ParaMCMC.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
