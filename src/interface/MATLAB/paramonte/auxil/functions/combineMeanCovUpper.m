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
%   we ask you to acknowledge the use of the ParaMonte library
%   in your work (education/research/industry/development/...)
%   by citing the ParaMonte library as described on this page:
%
%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This subroutine is the same as combineCovMean, with the IMPORTANT difference that
% only the upper triangles and diagonals of the input covariance matrices need to be given by the user: CovMatUpper, CovMatUpperB
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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