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
!>  This include file contains the implementation of procedures in [pm_distPareto](@ref pm_distPareto).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        !%%%%%%%%%%%%%%%%%%%%%%%%
#if     getParetoLogPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, alpha < 0._RKC, SK_"@getParetoLogPDFNF(): The condition `alpha < 0._RKC` must hold. alpha = "//getStr(alpha))
#if     ALD_ENABLED
        logPDFNF = log(-alpha) - alpha * logMinX
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getParetoLogPDFNF(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        logPDFNF = log(-alpha) - getLogSubExp(smaller = alpha * logMaxX, larger = alpha * logMinX)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getParetoLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, logMinX <= logx, SK_"@getParetoLogPDF(): The condition `logMinX <= logx` must hold. logx, logMaxX = "//getStr([logMinX, logx]))
#if     ALD_ENABLED
        call setParetoLogPDF(logPDF, logx, alpha, getParetoLogPDFNF(alpha, logMinX))
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logx <= logMaxX, SK_"@getParetoLogPDF(): The condition `logx <= logMaxX` must hold. logMinX, logx = "//getStr([logx, logMaxX]))
        call setParetoLogPDF(logPDF, logx, alpha, getParetoLogPDFNF(alpha, logMinX, logMaxX))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setParetoLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, alpha < 0._RKC, SK_"@setParetoLogPDF(): The condition `alpha < 0._RKC` must hold. alpha = "//getStr(alpha))
        logPDF = logPDFNF + (alpha - 1._RKC) * logx

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getParetoLogCDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

#if     ALD_ENABLED
        CHECK_ASSERTION(__LINE__, alpha < 0._RKC, SK_"@getParetoLogCDFNF(): The condition `alpha < 0._RKC` must hold. alpha = "//getStr(alpha))
        logCDFNF = 0._RKC
#elif   ALL_ENABLED
        real(RKC) :: alphaLogMinX
        CHECK_ASSERTION(__LINE__, alpha < 0._RKC, SK_"@getParetoLogCDFNF(): The condition `alpha < 0._RKC` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getParetoLogCDFNF(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        alphaLogMinX = alpha * logMinX
        logCDFNF = alphaLogMinX - getLogSubExp(smaller = alpha * logMaxX, larger = alphaLogMinX)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getParetoLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, logMinX <= logx, SK_"@getParetoLogCDF(): The condition `logx <= logMinX` must hold. logMinX, logx = "//getStr([logMinX, logx]))
#if     ALD_ENABLED
        call setParetoLogCDF(logCDF, logx, alpha, logMinX)
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logx <= logMaxX, SK_"@getParetoLogCDF(): The condition `logx <= logMaxX` must hold. logx, logMaxX = "//getStr([logx, logMaxX]))
        call setParetoLogCDF(logCDF, logx, alpha, logMinX, getParetoLogCDFNF(alpha, logMinX, logMaxX))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setParetoLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, alpha < 0._RKC, SK_"@setParetoLogCDF(): The condition `alpha < 0._RKC` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, logMinX <= logx, SK_"@setParetoLogCDF(): The condition `logMinX <= logx` must hold. logMinX, logx = "//getStr([logMinX, logx]))
        if (logx /= logMinX) then
            logCDF = getLogSubExp(smaller = alpha * (logx - logMinX), larger = 0._RKC)
#if         ALD_ENABLED
            CHECK_ASSERTION(__LINE__, exp(logCDF) <= 1._RKC + sqrt(epsilon(0._RKC)), \
            SK_"@setParetoLogCDF(): The condition `logCDF <= 0._RKC` must hold. The input arguments are inconsistent or `logx` is out of support. logCDF, logx, alpha, logMinX = "//getStr([logCDF, logx, alpha, logMinX]))
#elif       ALL_ENABLED
            logCDF = logCDF + logCDFNF
            CHECK_ASSERTION(__LINE__, exp(logCDF) <= 1._RKC + sqrt(epsilon(0._RKC)), \
            SK_"@setParetoLogCDF(): The condition `logCDF <= 0._RKC` must hold. The input arguments are inconsistent or `logx` is out of support. logCDF, logx, alpha, logMinX, logCDFNF = "//getStr([logCDF, logx, alpha, logMinX, logCDFNF]))
#else
#error      "Unrecognized interface."
#endif
        else
            logCDF = -log(huge(0._RKC))
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   getParetoLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

#if     ALD_ENABLED
        call setParetoLogQuan(logx, logCDF, alpha, logMinX)
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getParetoLogQuan(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        call setParetoLogQuan(logx, logCDF, alpha, logMinX, getParetoLogCDFNF(alpha, logMinX, logMaxX))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setParetoLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, alpha < 0._RKC, SK_"@setParetoLogQuan(): The condition `alpha < 0._RKC` must hold. alpha = "//getStr(alpha))
#if     LLALD_ENABLED
        CHECK_ASSERTION(__LINE__, logCDF < 0._RKC, SK_"@setParetoLogQuan(): The condition `logCDF < 0._RKC` must hold. logCDF = "//getStr(logCDF))
        logx = logMinX + getLogSubExp(smaller = logCDF, larger = 0._RKC) / alpha
#elif   LLALL_ENABLED
        CHECK_ASSERTION(__LINE__, logCDF <= 0._RKC, SK_"@setParetoLogQuan(): The condition `logCDF < 0._RKC` must hold. logCDF = "//getStr(logCDF))
        logx = logMinX + getLogSubExp(smaller = logCDF - logCDFNF, larger = 0._RKC) / alpha
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   getParetoLogRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

#if     ALD_ENABLED
        !call setParetoLogRand(logRand, getNegExpRand(sigma = 1._RKC), alpha, logMinX)
        call setParetoLogRand(logRand, log(1._RKC - getUnifRand(0._RKC, 1._RKC)), alpha, logMinX)
#elif   ALL_ENABLED
        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getParetoLogRand(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        !call setParetoLogRand(logRand, getNegExpRand(sigma = 1._RKC), alpha, logMinX, getParetoLogCDFNF(alpha, logMinX, logMaxX))
        call setParetoLogRand(logRand, log(1._RKC - getUnifRand(0._RKC, 1._RKC)), alpha, logMinX, getParetoLogCDFNF(alpha, logMinX, logMaxX))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setParetoLogRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, alpha < 0._RKC, SK_"@setParetoLogRand(): The condition `alpha < 0._RKC` must hold. alpha = "//getStr(alpha))
#if     LNALD_ENABLED
        CHECK_ASSERTION(__LINE__, negExpRand < 0._RKC, SK_"@setParetoLogRand(): The condition `negExpRand < 0._RKC` must hold. alpha = "//getStr(negExpRand))
        call setParetoLogQuan(logRand, negExpRand, alpha, logMinX)
#elif   LNALL_ENABLED
        CHECK_ASSERTION(__LINE__, negExpRand <= 0._RKC, SK_"@setParetoLogRand(): The condition `negExpRand <= 0._RKC` must hold. alpha = "//getStr(negExpRand))
        call setParetoLogQuan(logRand, negExpRand, alpha, logMinX, logCDFNF)
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif