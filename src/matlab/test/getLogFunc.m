%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                                                                            %%%%
%%%%    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           %%%%
%%%%                                                                                                                            %%%%
%%%%    Copyright (C) 2012-present, The Computational Data Science Lab                                                          %%%%
%%%%                                                                                                                            %%%%
%%%%    This file is part of the ParaMonte library.                                                                             %%%%
%%%%                                                                                                                            %%%%
%%%%    LICENSE                                                                                                                 %%%%
%%%%                                                                                                                            %%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          %%%%
%%%%                                                                                                                            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function logFunc = getLogFunc(point)

% This function returns the probability density function of the standard multivariate normal distribution of ndim dimensions.
    
    % Standard MultiVariate Normal (SMVN) specifications: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

%----------------------------------------------------------------------   

%   LOG_SMVN_COEF   = length(Point) * log( 1/sqrt(2*pi) );  
%   logFunc = LOG_SMVN_COEF - 0.5 * sum(Point.^2);

%----------------------------------------------------------------------

    logFunc = -0.5 * sum(point.^2);

%----------------------------------------------------------------------

    % Rosenbrock test function
%    logFunc = -( 100.0 * (point(2)-point(1)^2)^2 + (point(1)-1.0)^2 )/20;
%----------------------------------------------------------------------   

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  PYTHON-getLogFunc  %%%%%%%%%%%%%%%%%%%%%%
%%%%                                                    %%%%
%%%%  def getLogFunc(Point):                            %%%%
%%%%                                                    %%%%
%%%%      return -0.5 * np.sum( np.double( Point )**2 ) %%%%
%%%%                                                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
