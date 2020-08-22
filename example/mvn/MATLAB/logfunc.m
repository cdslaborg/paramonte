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
classdef logfunc

    properties (Constant)

        NDIM    = 4;                        % number of dimensions of the distribution
        MEAN    =   [ 0.0,0.0,0.0,0.0]';    % mean of the Multivariate Normal distribution
        COVMAT  =   [ 1.0,0.5,0.5,0.5       ... covariance matrix of the Multivariate Normal distribution
                    ; 0.5,1.0,0.5,0.5       ...
                    ; 0.5,0.5,1.0,0.5       ...
                    ; 0.5,0.5,0.5,1.0 ];
        INVCOV = inv(logfunc.COVMAT);       % the inverse of the covariance matrix

        % this is the log of the coefficient used in the definition of the MVN.

        MVN_COEF = logfunc.NDIM * log( 1. / sqrt(2.*pi) ) + log( sqrt(det(logfunc.INVCOV)) );

    end

    % define the objective function: getLogFunc

    methods (Static)
        function logFuncVal = get(point)
            % Return the natural logarithm of an NDIM-dimensional Multivariate Normal distribution
            % with the mean and covariance matrix as given in the above.
            % Reference: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
            %
            % note that the input point is an array of NDIM rows and 1 column.
            % Therefore, the MEAN array must have the same dimensions.

            % the logarithm of objective function: log(ObjectiveFunction)

            normedPoint = logfunc.MEAN - point;
            logFuncVal = logfunc.MVN_COEF - 0.5 * ( dot(normedPoint', logfunc.INVCOV * normedPoint) );

        end
    end

end
