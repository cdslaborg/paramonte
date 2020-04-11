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