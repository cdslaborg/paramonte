!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This include file contains the implementation of procedures in [pm_distPower](@ref pm_distPower).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%
#if     getPowerLogPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, alpha > 0._RKC, SK_"@getPowerLogPDFNF(): The condition `alpha > 0._RKC` must hold. alpha = "//getStr(alpha))
#if     ALD_ENABLED
        logPDFNF = log(alpha) - alpha * logMaxX
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getPowerLogPDFNF(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        logPDFNF = log(alpha) - getLogSubExp(smaller = alpha * logMinX, larger = alpha * logMaxX)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getPowerLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, logx <= logMaxX, SK_"@getPowerLogPDF(): The condition `logx <= logMaxX` must hold. logx, logMaxX = "//getStr([logx, logMaxX]))
#if     ALD_ENABLED
        call setPowerLogPDF(logPDF, logx, alpha, getPowerLogPDFNF(alpha, logMaxX))
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logMinX <= logx, SK_"@getPowerLogPDF(): The condition `logMinX <= logx` must hold. logMinX, logx = "//getStr([logMinX, logx]))
        call setPowerLogPDF(logPDF, logx, alpha, getPowerLogPDFNF(alpha, logMinX, logMaxX))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setPowerLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, alpha > 0._RKC, SK_"@setPowerLogPDF(): The condition `alpha > 0._RKC` must hold. alpha = "//getStr(alpha))
        logPDF = logPDFNF + (alpha - 1._RKC) * logx

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowerLogCDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

#if     ALD_ENABLED
        CHECK_ASSERTION(__LINE__, alpha > 0._RKC, SK_"@getPowerLogCDFNF(): The condition `alpha > 0._RKC` must hold. alpha = "//getStr(alpha))
        logCDFNF = - alpha * logMaxX
#elif   ALL_ENABLED
        real(RKC) :: alphaLogMinX
        alphaLogMinX = alpha * logMinX
        CHECK_ASSERTION(__LINE__, alpha > 0._RKC, SK_"@getPowerLogCDFNF(): The condition `alpha > 0._RKC` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getPowerLogCDFNF(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        logCDFNF = alphaLogMinX - getLogSubExp(smaller = alphaLogMinX, larger = alpha * logMaxX)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getPowerLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, logx <= logMaxX, SK_"@getPowerLogCDF(): The condition `logx <= logMaxX` must hold. logx, logMaxX = "//getStr([logx, logMaxX]))
#if     ALD_ENABLED
        call setPowerLogCDF(logCDF, logx, alpha, getPowerLogCDFNF(alpha, logMaxX))
#elif   ALL_ENABLED
        call setPowerLogCDF(logCDF, logx, alpha, logMinX, getPowerLogCDFNF(alpha, logMinX, logMaxX))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setPowerLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, alpha > 0._RKC, SK_"@setPowerLogCDF(): The condition `alpha > 0._RKC` must hold. alpha = "//getStr(alpha))
#if     ALD_ENABLED
        logCDF = logCDFNF + alpha * logx
        CHECK_ASSERTION(__LINE__, exp(logCDF) <= 1._RKC + sqrt(epsilon(0._RKC)), SK_"@setPowerLogCDF(): The condition `logCDF <= 0._RKC` must hold. The input arguments are inconsistent or `logx` is out of support. logCDF, logx, alpha, logCDFNF = "//getStr([logCDF, logx, alpha, logCDFNF]))
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logMinX <= logx, SK_"@setPowerLogCDF(): The condition `logMinX <= logx` must hold. logMinX, logx = "//getStr([logMinX, logx]))
        if (logx /= logMinX) then
            logCDF = logCDFNF + getLogSubExp(smaller = 0._RKC, larger = alpha * (logx - logMinX))
            CHECK_ASSERTION(__LINE__, exp(logCDF) <= 1._RKC + sqrt(epsilon(0._RKC)), SK_"@setPowerLogCDF(): The condition `logCDF <= 0._RKC` must hold. The input arguments are inconsistent or `logx` is out of support. epsilon(0._RKC), logCDF, logx, alpha, logMinX = "//getStr([epsilon(0._RKC), logCDF, logx, alpha, logMinX, logCDFNF]))
        else
            logCDF = -log(huge(0._RKC))
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowerLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

#if     ALD_ENABLED
        call setPowerLogQuan(logx, logCDF, alpha, getPowerLogCDFNF(alpha, logMaxX))
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getPowerLogQuan(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        call setPowerLogQuan(logx, logCDF, alpha, logMinX, getPowerLogCDFNF(alpha, logMinX, logMaxX))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowerLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, alpha > 0._RKC, SK_"@setPowerLogQuan(): The condition `alpha > 0._RKC` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, logCDF <= epsilon(0._RKC), SK_"@setPowerLogQuan(): The condition `logCDF <= epsilon(0._RKC)` must hold. logCDF, epsilon(0._RKC) = "//getStr([logCDF, epsilon(0._RKC)]))
#if         LLALD_ENABLED
            logx = (logCDF - logCDFNF) / alpha
#elif       LLALL_ENABLED
            logx = logMinX + getLogAddExp(getMinMax(logCDF - logCDFNF, 0._RKC)) / alpha
#else
#error      "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowerLogRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        use pm_distNegExp, only: getNegExpRand
        use pm_distPower, only: getPowerLogCDFNF
#if     ALD_ENABLED
        call setPowerLogRand(logRand, getNegExpRand(1._RKC), alpha, getPowerLogCDFNF(alpha, logMaxX))
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getPowerLogRand(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        call setPowerLogRand(logRand, getNegExpRand(1._RKC), alpha, logMinX, getPowerLogCDFNF(alpha, logMinX, logMaxX))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowerLogRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        use pm_distPower, only: setPowerLogQuan
        CHECK_ASSERTION(__LINE__, alpha > 0._RKC, SK_"@setPowerLogRand(): The condition `alpha > 0._RKC` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, negExpRand <= 0._RKC, SK_"@setPowerLogRand(): The condition `negExpRand <= 0._RKC` must hold. alpha = "//getStr(negExpRand))
#if     LNALD_ENABLED
        call setPowerLogQuan(logRand, negExpRand, alpha, logCDFNF)
#elif   LNALL_ENABLED
        call setPowerLogQuan(logRand, negExpRand, alpha, logMinX, logCDFNF)
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif