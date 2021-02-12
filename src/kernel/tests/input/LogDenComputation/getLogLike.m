function logLike = getLogLike(param)
    global DistSqSorted np LogFac
    nd = 2;
    logLam = param(1);
    logRad = param(2);
    logVol = log(pi) * nd * logRad;
    distSq = exp(2*logRad);
    logLPV = logLam + logVol;
    lamvol = exp(logLPV);
    LogProb = zeros(np,0);
    for i = 1:np
        count = length(find(DistSqSorted(:,i) <= distSq));
        LogProb(i) = count * logLPV - LogFac(count);
    end
    logLike = sum(LogProb) - np * lamvol;
end
