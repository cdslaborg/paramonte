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
% This subroutine is exact same thing as getSamCovMeanTrans, with the only differences that only the upper triangle of the covariance matrix is returned.
% Also, optional arguments are not available. This subroutine is specifically optimized for use in ParaMCMC.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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