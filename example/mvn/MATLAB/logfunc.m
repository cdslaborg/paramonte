%**********************************************************************************************************************************
%**********************************************************************************************************************************
%
%  ParaMonte: plain powerful parallel Monte Carlo library.
%
%  Copyright (C) 2012-present, The Computational Data Science Lab
%
%  This file is part of ParaMonte library. 
%
%  ParaMonte is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published by
%  the Free Software Foundation, version 3 of the License.
%
%  ParaMonte is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
%
%**********************************************************************************************************************************
%**********************************************************************************************************************************

classdef logfunc

    properties (Constant)
        NDIM    = 4; % number of dimensions of the distribution
        MEAN    =   [ 0.0,0.0,0.0,0.0]';    % mean of the Multivariate Normal distribution
        COVMAT  =   [ 1.0,0.5,0.5,0.5 ;     ... covariance matrix of the Multivariate Normal distribution
                    , 0.5,1.0,0.5,0.5 ;     ...
                    , 0.5,0.5,1.0,0.5 ;     ...
                    , 0.5,0.5,0.5,1.0 ];
    end

    % define the objective function: getLogFunc

    methods (Static)
        function logFunc = get(point)
            % note that the input point is an array of NDIM rows and 1 column.
            % Therefore, the MEAN array must have the same dimensions.
            logFunc = log(mvnpdf(point,logfunc.MEAN,logfunc.COVMAT));
        end
    end

end
