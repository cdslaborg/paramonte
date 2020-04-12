function burninLoc = getBurninLoc(lenLogFunc, refLogFunc, LogFunc)
    negLogIncidenceProb = log(lenLogFunc);
    burninLoc = 0;
    while true
        burninLoc = burninLoc + 1;
        if burninLoc < lenLogFunc && (refLogFunc - LogFunc(burninLoc)) > negLogIncidenceProb, continue; end
        break;
    end
end
