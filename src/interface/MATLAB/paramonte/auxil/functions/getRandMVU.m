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
function RandMVU = getRandMVU(nd, MeanVec, CholeskyLower)
    % Amir Shahmoradi, April 23, 2017, 1:36 AM, ICES, UTEXAS
    % Given the Cholesky Lower triangle and diagonals of a given covariance matrix, this function return one point
    % uniformly randomly drawn from inside of an nd-ellipsoid, whose nd elements getRandMVU(i), i=1,nd are guaranteed to be in range
    % MeanVec(i) - sqrt(CovMat(i,i)) < getRandMVU(i) < MeanVec(i) + sqrt(CovMat(i,i))

    DummyVec        = randn(nd,1);
    sumSqDummyVec   = sum(DummyVec.^2);
    
    dummy           = rand;
    dummy           = (dummy^(1/nd)) / sqrt(sumSqDummyVec);
    
    DummyVec        = DummyVec * dummy; % DummyVec(j) * dummy is a uniform random point from inside of nd-sphere
    
    RandMVU         = MeanVec' + CholeskyLower*DummyVec;
    RandMVU         = RandMVU';
    
end