function logDen = getLogDen(nd, logMeanMinDist, logVolUnitBall)
        logDen = nd * (gammaln((nd + 1) / nd) - logMeanMinDist) - logVolUnitBall;
end