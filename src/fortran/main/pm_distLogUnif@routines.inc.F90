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
!>  This include file contains the implementation of procedures in [pm_distLogUnif](@ref pm_distLogUnif).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%
#if     getLogUnifPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getLogUnifPDFNF(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        pdfnf = 1._RKG / (logMaxX - logMinX)

        !%%%%%%%%%%%%%%%%%%%%
#elif   getLogUnifPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKG < minx, SK_"@getLogUnifPDF(): The condition `0._RKG < minx` must hold. minx = "//getStr(minx))
        CHECK_ASSERTION(__LINE__, minx <= x, SK_"@getLogUnifPDF(): The condition `minx <= x` must hold. minx, x = "//getStr([minx, x]))
        CHECK_ASSERTION(__LINE__, x <= maxx, SK_"@getLogUnifPDF(): The condition `x <= maxx` must hold. x, exp(logMaxX) = "//getStr([x, maxx]))
        call setLogUnifPDF(pdf, x, getLogUnifPDFNF(log(minx), log(maxx)))

        !%%%%%%%%%%%%%%%%%%%%
#elif   setLogUnifPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKG <= x, SK_"@setLogUnifPDF(): The condition `0._RKG <= x` must hold. x = "//getStr(x))
        CHECK_ASSERTION(__LINE__, 0._RKG < pdfnf, SK_"@setLogUnifPDF(): The condition `0._RKG <= pdfnf` must hold. pdfnf = "//getStr(pdfnf))
        pdf = pdfnf / x

        !%%%%%%%%%%%%%%%%%%%%
#elif   getLogUnifCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, logMinX <= logx, SK_"@getLogUnifCDF(): The condition `logx <= logMinX` must hold. logMinX, logx = "//getStr([logMinX, logx]))
        CHECK_ASSERTION(__LINE__, logx <= logMaxX, SK_"@getLogUnifCDF(): The condition `logx <= logMaxX` must hold. logx, logMaxX = "//getStr([logx, logMaxX]))
        call setLogUnifCDF(cdf, logx, logMinX, getLogUnifPDFNF(logMinX, logMaxX))

        !%%%%%%%%%%%%%%%%%%%%
#elif   setLogUnifCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, logMinX <= logx, SK_"@setLogUnifCDF(): The condition `logMinX <= logx` must hold. logMinX, logx = "//getStr([logMinX, logx]))
        if (logMinX < logx) then
            cdf = (logx - logMinX) * pdfnf
            CHECK_ASSERTION(__LINE__, cdf < 1._RKG + sqrt(epsilon(0._RKG)), SK_"@setLogUnifCDF(): The condition `cdf <= 1._RKG` must hold. The input arguments are inconsistent or `logx` is out of support. cdf, logx, logMinX = "//getStr([cdf, logx, logMinX]))
        else
            cdf = 0._RKG
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getLogUnifLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, logMinX < logMaxX, SK_"@getLogUnifLogQuan(): The condition `logMinX < logMaxX` must hold. logMinX, logMaxX = "//getStr([logMinX, logMaxX]))
        call setLogUnifLogQuan(logx, cdf, logMinX, getLogUnifPDFNF(logMinX, logMaxX))

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setLogUnifLogQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKG <= cdf, SK_"@setLogUnifLogQuan(): The condition `0._RKG <= cdf` must hold. cdf = "//getStr(cdf))
        CHECK_ASSERTION(__LINE__, 1._RKG >= cdf, SK_"@setLogUnifLogQuan(): The condition `1._RKG >= cdf` must hold. cdf = "//getStr(cdf))
        logx = logMinX + cdf / pdfnf
        CHECK_ASSERTION(__LINE__, logMinX <= logx, SK_"@setLogUnifLogQuan(): The condition `logMinX <= logx` must hold. The input parameters are inconsistent. logMinX, logx = "//getStr([logMinX, logx]))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getLogUnifRand_ENABLED && MM_ENABLED && IK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer, parameter :: RKG = selected_real_kind(r = range(rand), radix = radix(rand))
        rand = int(getLogUnifRand(real(minx, RKG), real(maxx, RKG)), kind = IKG)
        !use pm_kind, only: RKS, RKD, RKQ, RKHR
        !if (real(huge(rand), RKH) < real(huge(0._RKS), RKH)) then
        !    rand = int(getLogUnifRand(real(minx, RKS), real(maxx, RKS)), kind = IKG)
        !    return
        !elseif (real(huge(rand), RKH) < real(huge(0._RKD), RKH)) then
        !    rand = int(getLogUnifRand(real(minx, RKD), real(maxx, RKD)), kind = IKG)
        !    return
        !else
        !    if (0 < RKH) then ! LCOV_EXCL_LINE
        !        if (real(huge(rand), RKH) < huge(0._RKH)) then ! LCOV_EXCL_LINE
        !            rand = int(getLogUnifRand(real(minx, RKH), real(maxx, RKH)), kind = IKG) ! LCOV_EXCL_LINE
        !            return ! LCOV_EXCL_LINE
        !        end if
        !    end if
        !end if
        !error stop "The kind of the input integer requires real numbers that can handle much larger values than what is supported by the available single, double, or quad precision real kinds. You must be living in the year 3000!" ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getLogUnifRand_ENABLED && MM_ENABLED && RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: logMinX, logRand
        CHECK_ASSERTION(__LINE__, minx < maxx, SK_"@getLogUnifRand(): The condition `minx < maxx` must hold. minx, maxx = "//getStr([minx, maxx]))
        logMinX = log(minx)
        call setLogUnifLogRand(logRand, getUnifRand(0._RKG, 1._RKG), logMinX, getLogUnifPDFNF(logMinX, log(maxx)))
        rand = exp(logRand)

        !%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setLogUnifLogRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKG <= urand, SK_"@setLogUnifLogRand(): The condition `0._RKG <= urand` must hold. urand = "//getStr(urand))
        CHECK_ASSERTION(__LINE__, 1._RKG >= urand, SK_"@setLogUnifLogRand(): The condition `1._RKG >= urand` must hold. urand = "//getStr(urand))
        logRand = logMinX + urand / pdfnf
        CHECK_ASSERTION(__LINE__, logMinX <= logRand, SK_"@setLogUnifLogRand(): The condition `logMinX <= logRand` must hold. The input parameters are inconsistent. logMinX, logRand = "//getStr([logMinX, logRand]))

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif