% returns log( exp(logFunc1) - exp(logFunc2) ) robustly
% onus is on the user to ensure logFunc1 > logFunc2. To ensure wrong input does not go unnoticed
function logSubExp = getLogSubExp(logFunc1, logFunc2)

    logSubExp = logFunc1 + log( 1 - exp(logFunc2 - logFunc1) );

end