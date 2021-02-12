function logDis = getLogDis(nd, logDen, logVolUnitBall)
        logDis = gammaln((nd + 1) / nd) - (logDen + logVolUnitBall) / nd;
end