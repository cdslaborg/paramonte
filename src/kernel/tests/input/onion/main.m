%clear all
%close all

logpi = log(pi);

logLambdaMin = log(1);
logLambdaMax = log(10000);
logLambdaDel = (logLambdaMax - logLambdaMin) / 19;
LogLam = logLambdaMin : logLambdaDel : logLambdaMax;
lenLogLambda = length(LogLam);
LogVolTrue = zeros(lenLogLambda,1);

ndRange = 128; [1,2,4,8,16,32,64,128];
ndRangeLen = length(ndRange);

LogVol = zeros(lenLogLambda, ndRangeLen);
for id = 1:ndRangeLen

    nd = ndRange(id)
    ndHalf = nd / 2.;
    logVolUnitBall = ndHalf * logpi - gammaln(ndHalf + 1);

    for i = 1:lenLogLambda
        lambda = exp(LogLam(i));
        np = max([3, poissrnd(lambda)]);
        Point = unifrnd(0, 1, np, nd);
        Dist = squareform(pdist(Point));
        DistSorted = sort(Dist, "ascend");
        logMeanMinDist = log(mean(DistSorted(2,:)));
        LogVol(i,id) = log(np) - getLogDen(nd, logMeanMinDist, logVolUnitBall);
    end

end
% figure; hold on;
% plot( Point(:,1) ...
%     , Point(:,2) ...
%     , '.' ...
%     );
% hold off;

for id = 1:ndRangeLen
    figure; hold on;
    plot( LogLam ...
        , LogVol(:,id) ...
        , '.' ...
        , "markersize", 10 ...
        );
    plot( LogLam ...
        , LogVolTrue ...
        , "-" ...
        , "linewidth", 3 ...
        );
    hold off;
end

%figure; hold on;
%histogram(LogDen - LogLam);
%hold off;

SumMinDist = zeros(np,1);
LogVolCorrected = zeros(np,1);
LogMeanMinDist = zeros(np,1);
SumMinDist(1) = sum(DistSorted(2,:));
SumDistSorted = sum(DistSorted);
npUppLim = np-2;
for ip = 1:npUppLim
    npRemained = np - ip;
    [maxval, maxloc] = max(SumDistSorted);
    SumDistSorted(maxloc(1)) = -inf;
    SumMinDist(ip+1) = SumMinDist(ip) - DistSorted(2,maxloc(1));
    LogMeanMinDist(ip) = log(SumMinDist(ip+1) / npRemained);
    LogVolCorrected(ip) = log(np) - getLogDen(nd, LogMeanMinDist(ip), logVolUnitBall);
end

figure; hold on;
plot( npUppLim:-1:1 ...
    , LogVolCorrected(1:npUppLim) ...
    , '.' ...
    , "markersize", 10 ...
    );
plot( npUppLim:-1:1 ...
    , zeros(npUppLim,1) ...
    , '.' ...
    , "markersize", 10 ...
    );
set(gca, "xscale", "log", "yscale", "linear")
hold off;

