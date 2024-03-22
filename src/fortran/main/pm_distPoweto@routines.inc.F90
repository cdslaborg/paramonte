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
!>  This include file contains the implementation of procedures in [pm_distPoweto](@ref pm_distPoweto).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%
#if     getPowetoLogPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0 < alpha .or. present(logMinX), SK_"@getPowetoLogPDFNF(): The condition `0 < alpha .or. present(logMinX)` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, alpha < 0 .or. present(logMaxX), SK_"@getPowetoLogPDFNF(): The condition `alpha < 0 .or. present(logMaxX)` must hold. alpha = "//getStr(alpha))

        if (alpha < 0._RKC) then
            if (present(logMaxX)) then
                logPDFNF = getParetoLogPDFNF(alpha, logMinX, logMaxX)
            else
                logPDFNF = getParetoLogPDFNF(alpha, logMinX)
            end if
        elseif (0._RKC < alpha) then
            if (present(logMinX)) then
                logPDFNF = getPowerLogPDFNF(alpha, logMinX, logMaxX)
            else
                logPDFNF = getPowerLogPDFNF(alpha, logMaxX)
            end if
        else
            logPDFNF = -log(logMaxX - logMinX)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowetoLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        call setPowetoLogPDF(logPDF, logx, alpha, getPowetoLogPDFNF(alpha, logMinX, logMaxX))

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowetoLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        if (0._RKC < alpha) then
            call setParetoLogPDF(logPDF, logx, alpha, logPDFNF)
        elseif (alpha < 0._RKC) then
            call setPowerLogPDF(logPDF, logx, alpha, logPDFNF)
        else
            logPDF = logPDFNF - logx
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowetoLogCDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0 < alpha .or. present(logMinX), SK_"@getPowetoLogCDFNF(): The condition `0 < alpha .or. present(logMinX)` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, alpha < 0 .or. present(logMaxX), SK_"@getPowetoLogCDFNF(): The condition `alpha < 0 .or. present(logMaxX)` must hold. alpha = "//getStr(alpha))

        if (alpha < 0._RKC) then
            if (present(logMaxX)) then
                logCDFNF = getParetoLogCDFNF(alpha, logMinX, logMaxX)
            else
                logCDFNF = getParetoLogCDFNF(alpha, logMinX)
            end if
        elseif (0._RKC < alpha) then
            if (present(logMinX)) then
                logCDFNF = getPowerLogCDFNF(alpha, logMinX, logMaxX)
            else
                logCDFNF = getPowerLogCDFNF(alpha, logMaxX)
            end if
        else
            logCDFNF = getPowetoLogPDFNF(alpha, logMinX, logMaxX)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowetoLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        call setPowetoLogCDF(logCDF, logx, alpha, logMinX, getPowetoLogCDFNF(alpha, logMinX, logMaxX))

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowetoLogCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0 < alpha .or. present(logMinX), SK_"@getPowetoLogCDFNF(): The condition `0 < alpha .or. present(logMinX)` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, alpha < 0 .or. present(logCDFNF), SK_"@getPowetoLogCDFNF(): The condition `alpha < 0 .or. present(logCDFNF)` must hold. alpha = "//getStr(alpha))

        if (alpha < 0._RKC) then
            if (present(logCDFNF)) then
                call setParetoLogCDF(logCDF, logx, alpha, logMinX, logCDFNF)
            else
                call setParetoLogCDF(logCDF, logx, alpha, logMinX)
            end if
        elseif (0._RKC < alpha) then
            if (present(logMinX)) then
                call setPowerLogCDF(logCDF, logx, alpha, logMinX, logCDFNF)
            else
                call setPowerLogCDF(logCDF, logx, alpha, logCDFNF)
            end if
        else
            logCDF = log(logx - logMinX) + logCDFNF
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowetoLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        call setPowetoLogQuan(logx, logCDF, alpha, logMinX, getPowetoLogCDFNF(alpha, logMinX, logMaxX))

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowetoLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0 < alpha .or. present(logMinX), SK_"@getPowetoLogCDFNF(): The condition `0 < alpha .or. present(logMinX)` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, alpha < 0 .or. present(logCDFNF), SK_"@getPowetoLogCDFNF(): The condition `alpha < 0 .or. present(logCDFNF)` must hold. alpha = "//getStr(alpha))

        if (alpha < 0._RKC) then
            if (present(logCDFNF)) then
                call setParetoLogQuan(logx, logCDF, alpha, logMinX, logCDFNF)
            else
                call setParetoLogQuan(logx, logCDF, alpha, logMinX)
            end if
        elseif (0._RKC < alpha) then
            if (present(logMinX)) then
                call setPowerLogQuan(logx, logCDF, alpha, logMinX, logCDFNF)
            else
                call setPowerLogQuan(logx, logCDF, alpha, logCDFNF)
            end if
        else
            logx = exp(logCDF - logCDFNF) + logMinX
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowetoLogRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        call setPowetoLogRand(logRand, getNegExpRand(1._RKC), alpha, logMinX, getPowetoLogCDFNF(alpha, logMinX, logMaxX))

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowetoLogRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, negExpRand <= 0._RKC, SK_"@setPowetoLogRand(): The condition `negExpRand <= 0._RKC` must hold. alpha = "//getStr(negExpRand))
        call setPowetoLogQuan(logRand, negExpRand, alpha, logMinX, logCDFNF)

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif