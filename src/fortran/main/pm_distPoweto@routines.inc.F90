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
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%
#if     getPowetoLogPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%
        use pm_mathCumSum, only: setCumSum
        use pm_mathLogSubExp, only: getLogSubExp
        use pm_mathLogSumExp, only: getLogSumExp
        real(RKC), parameter :: LOG_HUGE = log(huge(LogLimX))
        integer(IK) :: i, lenAlpha, lenLogLimX
        real(RKC)   :: maxArea, logIntegral
#if     ALD_ENABLED
        real(RKC)   :: CumSumArea(size(LogLimX, kind = IK)) ! initially, LogArea.
#elif   ALC_ENABLED
        CHECK_ASSERTION(__LINE__, size(LogLimX,1,IK) == size(CumSumArea,1,IK), SK_"@getPowetoLogPDFNF(): The condition `size(LogLimX,1,IK) == size(CumSumArea,1,IK)` must hold. size(LogLimX), size(CumSumArea) = "//getStr([size(LogLimX,1,IK), size(CumSumArea,1,IK)]))
#else
#error  "Unrecognized interface."
#endif
        lenAlpha = size(Alpha, kind = IK)
        lenLogLimX = size(LogLimX, kind = IK)
        CHECK_ASSERTION(__LINE__, lenAlpha > 0_IK, SK_"@getPowetoLogPDFNF(): The condition `size(Alpha) > 0` must hold. size(Alpha) = "//getStr(lenAlpha))
        CHECK_ASSERTION(__LINE__, lenLogLimX == lenAlpha + 1_IK, SK_"@getPowetoLogPDFNF(): The condition `size(LogLimX) == size(Alpha) + 1` must hold. size(LogLimX), size(Alpha) = "//getStr([lenLogLimX, lenAlpha]))
        CHECK_ASSERTION(__LINE__, isAscending(LogLimX), SK_"@getPowetoLogPDFNF(): The condition `isAscending(LogLimX)` must hold. LogLimX = "//getStr(LogLimX))
        CHECK_ASSERTION(__LINE__, Alpha(1) > 0._RKC .or. LogLimX(1) > -Log_HUGE, SK_"@getPowetoLogPDFNF(): The conditions `Alpha(1) > 0._RKC .or. LogLimX(1) > -Log_HUGE` must hold. Alpha(1), LogLimX(1), -Log_HUGE = "//getStr([Alpha(1), LogLimX(1), -Log_HUGE]))
        CHECK_ASSERTION(__LINE__, Alpha(lenAlpha) < 0._RKC .or. LogLimX(lenLogLimX) < Log_HUGE, SK_"@getPowetoLogPDFNF(): The conditions `Alpha(lenAlpha) < 0._RKC .or. LogLimX(lenLogLimX) < Log_HUGE` must hold. Alpha(lenAlpha), LogLimX(lenLogLimX), Log_HUGE = "//getStr([Alpha(lenAlpha), LogLimX(lenLogLimX), Log_HUGE]))
        ! Initially, CumSumArea contains areas of the segments.
        if (Alpha(lenAlpha) > 0._RKC) then ! Power tail.
            CumSumArea(lenLogLimX) = getLogSubExp(smaller = Alpha(lenAlpha) * LogLimX(lenAlpha), larger = Alpha(lenAlpha) * LogLimX(lenLogLimX)) - log(Alpha(lenAlpha))
        elseif (Alpha(lenAlpha) < 0._RKC) then
            CumSumArea(lenLogLimX) = getLogSubExp(smaller = Alpha(lenAlpha) * LogLimX(lenLogLimX), larger = Alpha(lenAlpha) * LogLimX(lenAlpha)) - log(-Alpha(lenAlpha))
        else
            CumSumArea(lenLogLimX) = log(LogLimX(lenLogLimX) - LogLimX(lenAlpha))
        end if
        LogNormFac(lenAlpha) = 0._RKC
        maxArea = CumSumArea(lenLogLimX)
        do i = lenAlpha - 1, 1_IK, -1_IK
            LogNormFac(i) = LogNormFac(i + 1) + (Alpha(i + 1) - Alpha(i)) * LogLimX(i + 1)
            if (Alpha(i) > 0._RKC) then
                CumSumArea(i + 1) = LogNormFac(i) + getLogSubExp(smaller = Alpha(i) * LogLimX(i), larger = Alpha(i) * LogLimX(i + 1)) - log(Alpha(i))
            elseif (Alpha(i) < 0._RKC) then
                CumSumArea(i + 1) = LogNormFac(i) + getLogSubExp(smaller = Alpha(i) * LogLimX(i + 1), larger = Alpha(i) * LogLimX(i)) - log(-Alpha(i))
            else
                CumSumArea(i + 1) = LogNormFac(i) + log(LogLimX(i + 1) - LogLimX(i))
            end if
            if (maxArea < CumSumArea(i + 1)) maxArea = CumSumArea(i + 1)
        end do
        logIntegral = getLogSumExp(CumSumArea(2:lenLogLimX), maxArea)
        LogNormFac = LogNormFac - logIntegral
#if     ALC_ENABLED
        CumSumArea(1) = 0._RKC
        CumSumArea(2:lenLogLimX) = exp(CumSumArea(2:lenLogLimX) - logIntegral)
        call setCumSum(CumSumArea(2:lenLogLimX))
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowetoLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        if (present(LogNormFac)) then
            call setPowetoLogPDF(logPDF, logx, Alpha, LogLimX, LogNormFac)
        else
            call setPowetoLogPDF(logPDF, logx, Alpha, LogLimX, getPowetoLogPDFNF(Alpha, LogLimX))
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowetoLogPDF_ENABLED && ALL_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenAlpha
        lenAlpha = size(Alpha, kind = IK)
        CHECK_ASSERTION(__LINE__, lenAlpha > 0_IK, SK_"@setPowetoLogPDF(): The condition `size(Alpha) > 0` must hold. size(LogLimX) = "//getStr(lenAlpha))
        CHECK_ASSERTION(__LINE__, size(LogNormFac,1,IK) == lenAlpha, SK_"@setPowetoLogPDF(): The condition `size(LogNormFac,1,IK) == lenAlpha` must hold. size(LogNormFac,1,IK), lenAlpha = "//getStr([size(LogNormFac,1,IK), lenAlpha]))
        CHECK_ASSERTION(__LINE__, size(LogLimX,1,IK) == lenAlpha + 1_IK, SK_"@setPowetoLogPDF(): The condition `size(LogLimX) == size(Alpha) + 1` must hold. size(LogLimX,1,IK), size(Alpha,1,IK) = "//getStr([size(LogLimX,1,IK), size(Alpha,1,IK)]))
        CHECK_ASSERTION(__LINE__, isAscending(LogLimX), SK_"@setPowetoLogPDF(): The conditions `isAscending(LogLimX)` must hold. LogLimX = "//getStr(LogLimX))
        CHECK_ASSERTION(__LINE__, LogLimX(1) <= logx, SK_"@setPowetoLogPDFNF(): The condition `LogLimX(1) <= logx` must hold. LogLimX(1), logx = "//getStr([LogLimX(1), logx]))
        CHECK_ASSERTION(__LINE__, logx <= LogLimX(size(LogLimX,1,IK)), SK_"@getPowetoLogPDFNF(): The condition `logx < LogLimX(size(LogLimX,1,IK))` must hold. logx, LogLimX(size(LogLimX,1,IK)) = "//getStr([logx, LogLimX(size(LogLimX,1,IK))]))
        !CHECK_ASSERTION(__LINE__, LogLimX(size(LogLimX,1,IK)) < log(huge(LogLimX)), SK_"@getPowetoLogPDFNF(): The condition `LogLimX(size(LogLimX)) < log(huge(LogLimX))` must hold. LogLimX = "//getStr(LogLimX))
        do i = lenAlpha, 1_IK, -1_IK
            if (logx < LogLimX(i)) cycle
            logPDF = LogNormFac(i) + (Alpha(i) - 1._RKC) * logx
            return
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowetoLogPDF_ENABLED && BAN_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0_IK < bin .and. bin <= size(Alpha, kind = IK), SK_"@setPowetoLogPDFNF(): The conditions `0_IK < bin .and. bin <= size(Alpha, kind = IK)` must hold. bin, size(Alpha, kind = IK) = "//getStr([bin, size(Alpha, kind = IK)]))
        CHECK_ASSERTION(__LINE__, size(Alpha, kind = IK) == size(LogNormFac, kind = IK), SK_"@setPowetoLogPDFNF(): The condition `size(Alpha) == size(LogNormFac)` must hold. size(Alpha), size(LogNormFac) = "//getStr([size(Alpha, kind = IK), size(LogNormFac, kind = IK)]))
        logPDF = LogNormFac(bin) + (Alpha(bin) - 1._RKC) * logx

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowetoCDF_ENABLED && ALDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKC) :: CumSumArea(size(LogLimX, kind = IK))
        call setPowetoCDF(cdf, logx, Alpha, LogLimX, getPowetoLogPDFNF(Alpha, LogLimX, CumSumArea), CumSumArea)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPowetoCDF_ENABLED && ALLC_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setPowetoCDF(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowetoCDF_ENABLED && BAN_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_mathLogSubExp, only: getLogSubExp
        CHECK_ASSERTION(__LINE__, size(LogLimX, kind = IK) == size(CumSumArea, kind = IK), SK_"@setPowetoCDF(): The condition `size(LogLimX) == size(CumSumArea)` must hold. size(LogLimX), size(CumSumArea) = "//getStr([size(LogLimX, kind = IK), size(CumSumArea, kind = IK)]))
        CHECK_ASSERTION(__LINE__, size(LogLimX, kind = IK) == size(Alpha, kind = IK) + 1_IK, SK_"@setPowetoCDF(): The condition `size(LogLimX) == size(Alpha) + 1` must hold. size(LogLimX), size(Alpha) = "//getStr([size(LogLimX, kind = IK), size(Alpha, kind = IK)]))
        CHECK_ASSERTION(__LINE__, size(Alpha, kind = IK) == size(LogNormFac, kind = IK), SK_"@setPowetoCDF(): The condition `size(Alpha) == size(LogNormFac)` must hold. size(Alpha), size(LogNormFac) = "//getStr([size(Alpha, kind = IK), size(LogNormFac, kind = IK)]))
        CHECK_ASSERTION(__LINE__, 0_IK < bin .and. bin < size(LogLimX, kind = IK), SK_"@setPowetoCDF(): The conditions `0_IK < bin .and. bin < size(LogLimX)` must hold. bin, size(LogLimX) = "//getStr([bin, size(LogLimX, kind = IK)]))
        CHECK_ASSERTION(__LINE__, LogLimX(bin) <= logx .and. logx <= LogLimX(bin + 1), SK_"@setPowetoCDF(): The conditions `LogLimX(bin) <= logx .and. logx < LogLimX(bin + 1)` must hold. LogLimX(bin), logx, LogLimX(bin+1) = "//getStr([LogLimX(bin), logx, LogLimX(bin+1)]))
#if     CHECK_ENABLED
        block
            real(RKC) :: LogNormFacRef(size(LogNormFac, kind = IK)), CumSumAreaRef(size(CumSumArea, kind = IK))
            LogNormFacRef = getPowetoLogPDFNF(Alpha, LogLimX, CumSumAreaRef)
            CHECK_ASSERTION(__LINE__, all(abs(LogNormFac - LogNormFacRef) <= epsilon(0._RKC) * 1000), SK_"@setPowetoCDF(): The condition `all(abs(LogNormFac - LogNormFacRef) <= epsilon(0._RKC) * 1000)` must hold. LogNormFac - LogNormFacRef = "//getStr([LogNormFac - LogNormFacRef]))
            CHECK_ASSERTION(__LINE__, all(abs(CumSumArea - CumSumAreaRef) <= epsilon(0._RKC) * 1000), SK_"@setPowetoCDF(): The condition `all(abs(CumSumArea - CumSumAreaRef) <= epsilon(0._RKC) * 1000)` must hold. CumSumArea - CumSumAreaRef = "//getStr([CumSumArea - CumSumAreaRef]))
        end block
#endif
        cdf = CumSumArea(bin)
        if (logx /= LogLimX(bin)) then
            if (Alpha(bin) > 0._RKC) then
                cdf = CumSumArea(bin) + exp(LogNormFac(bin) + getLogSubExp(smaller = Alpha(bin) * LogLimX(bin), larger = Alpha(bin) * logx)) / Alpha(bin)
            elseif (Alpha(bin) < 0._RKC) then
                cdf = CumSumArea(bin) - exp(LogNormFac(bin) + getLogSubExp(smaller = Alpha(bin) * logx, larger = Alpha(bin) * LogLimX(bin))) / Alpha(bin)
            else
                cdf = CumSumArea(bin) + exp(LogNormFac(bin) + log(logx - LogLimX(bin)))
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPowetoCDF_ENABLED && MAN_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenAlpha
        lenAlpha = size(Alpha, kind = IK)
        CHECK_ASSERTION(__LINE__, lenAlpha > 0_IK, SK_"@setPowetoCDF(): The condition `size(Alpha) > 0` must hold. size(LogLimX) = "//getStr(lenAlpha))
        CHECK_ASSERTION(__LINE__, size(LogNormFac,1,IK) == lenAlpha, SK_"@setPowetoCDF(): The condition `size(LogNormFac,1,IK) == lenAlpha` must hold. size(LogNormFac,1,IK), lenAlpha = "//getStr([size(LogNormFac,1,IK), lenAlpha]))
        CHECK_ASSERTION(__LINE__, size(LogLimX,1,IK) == lenAlpha + 1_IK, SK_"@setPowetoCDF(): The condition `size(LogLimX) == size(Alpha) + 1` must hold. size(LogLimX,1,IK), size(Alpha,1,IK) = "//getStr([size(LogLimX,1,IK), size(Alpha,1,IK)]))
        CHECK_ASSERTION(__LINE__, LogLimX(1) <= logx .and. logx <= LogLimX(lenAlpha + 1), SK_"@setPowetoCDF(): The condition `LogLimX(1) <= logx and. logx < LogLimX(size(LogLimX))` must hold. LogLimX(1), logx, LogLimX(size(LogLimX)) = "//getStr([LogLimX(1), logx, LogLimX(lenAlpha + 1)]))
        CHECK_ASSERTION(__LINE__, isAscending(LogLimX), SK_"@setPowetoCDF(): The conditions `isAscending(LogLimX)` must hold. LogLimX = "//getStr(LogLimX))
        !CHECK_ASSERTION(__LINE__, LogLimX(size(LogLimX,1,IK)) < log(huge(LogLimX)), SK_"@getPowetoLogPDFNF(): The condition `LogLimX(size(LogLimX)) < log(huge(LogLimX))` must hold. LogLimX = "//getStr(LogLimX))
        !CHECK_ASSERTION(__LINE__, merge(Alpha(lenAlpha) < -1._RKC, .true., size(LogLimX,1,IK) == lenAlpha), SK_"@getPowetoLogPDFNF(): The condition `Alpha(size(Alpha)) < -1._RKC` must hold. Alpha = "//getStr(Alpha))
        do i = lenAlpha, 1_IK, -1_IK
            if (logx < LogLimX(i)) cycle
            call setPowetoCDF(cdf, logx, Alpha, LogLimX, LogNormFac, CumSumArea, i)
            return
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif