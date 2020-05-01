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
    %logFunc = -log( 100.0 * (point(2)-point(1)^2)^2 + (point(1)-1.0)^2 );
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
