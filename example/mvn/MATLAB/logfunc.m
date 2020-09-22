%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Description:
%%%%        +   Return the natural logarithm of an ndim-dimensional Multivariate Normal (MVN)
%%%%            probability density function (PDF) with the Mean and Covariance Matrix as defined below.
%%%%            Reference: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
%%%%    Input:
%%%%        +   point:      The input 64-bit real-valued vector of length ndim,
%%%%                        at which the natural logarithm of objective function is computed.
%%%%  Output:
%%%%        +   logFuncVal: A 64-bit real scalar number representing the natural logarithm of the objective function.
%%%%  Author:
%%%%        +   Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
%%%%  Visit:
%%%%        +   https://www.cdslab.org/paramonte
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
