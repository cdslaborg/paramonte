function RandMVN = getRandMVN(nd, MeanVec, CholeskyLower)
    % Amir Shahmoradi, April 23, 2017, 12:36 AM, ICES, UTEXAS
    % returns a multivariate Normal random vector.

    % CholeskyLower, Diagonal: Cholesky lower triangle and its diagonal terms, calculated from the input CovMat.
    
    RandMVN = MeanVec' + CholeskyLower*randn(nd,1);
    RandMVN = RandMVN';

end