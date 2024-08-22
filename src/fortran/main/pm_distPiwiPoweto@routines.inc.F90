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
!>  This include file contains the implementation of procedures in [pm_distPiwiPoweto](@ref pm_distPiwiPoweto).
!>  
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getPiwiPowetoLogPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_mathCumSum, only: setCumSum
        use pm_mathLogSubExp, only: getLogSubExp
        use pm_mathLogSumExp, only: getLogSumExp
        real(RKG), parameter :: LOG_HUGE = log(huge(logLimX))
        integer(IK) :: i, lenAlpha, lenLogLimX
        real(RKG)   :: maxArea, logIntegral
#if     ALD_ENABLED
        real(RKG)   :: cumSumArea(size(logLimX, kind = IK)) ! initially, LogArea.
#elif   ALC_ENABLED
        CHECK_ASSERTION(__LINE__, size(logLimX,1,IK) == size(cumSumArea,1,IK), SK_"@getPiwiPowetoLogPDFNF(): The condition `size(logLimX,1,IK) == size(cumSumArea,1,IK)` must hold. size(logLimX), size(cumSumArea) = "//getStr([size(logLimX,1,IK), size(cumSumArea,1,IK)]))
#else
#error  "Unrecognized interface."
#endif
        lenAlpha = size(alpha, kind = IK)
        lenLogLimX = size(logLimX, kind = IK)
        CHECK_ASSERTION(__LINE__, lenAlpha > 0_IK, SK_"@getPiwiPowetoLogPDFNF(): The condition `size(alpha) > 0` must hold. size(alpha) = "//getStr(lenAlpha))
        CHECK_ASSERTION(__LINE__, lenLogLimX == lenAlpha + 1_IK, SK_"@getPiwiPowetoLogPDFNF(): The condition `size(logLimX) == size(alpha) + 1` must hold. size(logLimX), size(alpha) = "//getStr([lenLogLimX, lenAlpha]))
        CHECK_ASSERTION(__LINE__, isAscending(logLimX), SK_"@getPiwiPowetoLogPDFNF(): The condition `isAscending(logLimX)` must hold. logLimX = "//getStr(logLimX))
        CHECK_ASSERTION(__LINE__, alpha(1) > 0._RKG .or. logLimX(1) > -Log_HUGE, SK_"@getPiwiPowetoLogPDFNF(): The conditions `alpha(1) > 0._RKG .or. logLimX(1) > -Log_HUGE` must hold. alpha(1), logLimX(1), -Log_HUGE = "//getStr([alpha(1), logLimX(1), -Log_HUGE]))
        CHECK_ASSERTION(__LINE__, alpha(lenAlpha) < 0._RKG .or. logLimX(lenLogLimX) < Log_HUGE, SK_"@getPiwiPowetoLogPDFNF(): The conditions `alpha(lenAlpha) < 0._RKG .or. logLimX(lenLogLimX) < Log_HUGE` must hold. alpha(lenAlpha), logLimX(lenLogLimX), Log_HUGE = "//getStr([alpha(lenAlpha), logLimX(lenLogLimX), Log_HUGE]))

        ! Initially, cumSumArea contains areas of the segments.
        if (alpha(lenAlpha) > 0._RKG) then ! Power tail.
            cumSumArea(lenLogLimX) = getLogSubExp(smaller = alpha(lenAlpha) * logLimX(lenAlpha), larger = alpha(lenAlpha) * logLimX(lenLogLimX)) - log(alpha(lenAlpha))
        elseif (alpha(lenAlpha) < 0._RKG) then
            cumSumArea(lenLogLimX) = getLogSubExp(smaller = alpha(lenAlpha) * logLimX(lenLogLimX), larger = alpha(lenAlpha) * logLimX(lenAlpha)) - log(-alpha(lenAlpha))
        else
            cumSumArea(lenLogLimX) = log(logLimX(lenLogLimX) - logLimX(lenAlpha))
        end if
        logPDFNF(lenAlpha) = 0._RKG
        maxArea = cumSumArea(lenLogLimX)
        do i = lenAlpha - 1, 1_IK, -1_IK
            logPDFNF(i) = logPDFNF(i + 1) + (alpha(i + 1) - alpha(i)) * logLimX(i + 1)
            if (alpha(i) > 0._RKG) then
                cumSumArea(i + 1) = logPDFNF(i) + getLogSubExp(smaller = alpha(i) * logLimX(i), larger = alpha(i) * logLimX(i + 1)) - log(alpha(i))
            elseif (alpha(i) < 0._RKG) then
                cumSumArea(i + 1) = logPDFNF(i) + getLogSubExp(smaller = alpha(i) * logLimX(i + 1), larger = alpha(i) * logLimX(i)) - log(-alpha(i))
            else
                cumSumArea(i + 1) = logPDFNF(i) + log(logLimX(i + 1) - logLimX(i))
            end if
            if (maxArea < cumSumArea(i + 1)) maxArea = cumSumArea(i + 1)
        end do
        logIntegral = getLogSumExp(cumSumArea(2:lenLogLimX), maxArea)
        logPDFNF = logPDFNF - logIntegral
#if     ALC_ENABLED
        cumSumArea(1) = 0._RKG
        cumSumArea(2:lenLogLimX) = exp(cumSumArea(2:lenLogLimX) - logIntegral)
        call setCumSum(cumSumArea(2:lenLogLimX))
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPiwiPowetoLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(logPDFNF)) then
            call setPiwiPowetoLogPDF(logPDF, logx, alpha, logLimX, logPDFNF)
        else
            call setPiwiPowetoLogPDF(logPDF, logx, alpha, logLimX, getPiwiPowetoLogPDFNF(alpha, logLimX))
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPiwiPowetoLogPDF_ENABLED && ALL_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenAlpha
        lenAlpha = size(alpha, kind = IK)
        CHECK_ASSERTION(__LINE__, lenAlpha > 0_IK, SK_"@setPiwiPowetoLogPDF(): The condition `size(alpha) > 0` must hold. size(logLimX) = "//getStr(lenAlpha))
        CHECK_ASSERTION(__LINE__, size(logPDFNF,1,IK) == lenAlpha, SK_"@setPiwiPowetoLogPDF(): The condition `size(logPDFNF,1,IK) == lenAlpha` must hold. size(logPDFNF,1,IK), lenAlpha = "//getStr([size(logPDFNF,1,IK), lenAlpha]))
        CHECK_ASSERTION(__LINE__, size(logLimX,1,IK) == lenAlpha + 1_IK, SK_"@setPiwiPowetoLogPDF(): The condition `size(logLimX) == size(alpha) + 1` must hold. size(logLimX,1,IK), size(alpha,1,IK) = "//getStr([size(logLimX,1,IK), size(alpha,1,IK)]))
        CHECK_ASSERTION(__LINE__, isAscending(logLimX), SK_"@setPiwiPowetoLogPDF(): The conditions `isAscending(logLimX)` must hold. logLimX = "//getStr(logLimX))
        CHECK_ASSERTION(__LINE__, logLimX(1) <= logx, SK_"@setPiwiPowetoLogPDFNF(): The condition `logLimX(1) <= logx` must hold. logLimX(1), logx = "//getStr([logLimX(1), logx]))
        CHECK_ASSERTION(__LINE__, logx <= logLimX(size(logLimX,1,IK)), SK_"@getPiwiPowetoLogPDFNF(): The condition `logx < logLimX(size(logLimX,1,IK))` must hold. logx, logLimX(size(logLimX,1,IK)) = "//getStr([logx, logLimX(size(logLimX,1,IK))]))
        !check_assertion(__LINE__, logLimX(size(logLimX,1,IK)) < log(huge(logLimX)), SK_"@getPiwiPowetoLogPDFNF(): The condition `logLimX(size(logLimX)) < log(huge(logLimX))` must hold. logLimX = "//getStr(logLimX))
        do i = lenAlpha, 1_IK, -1_IK
            if (logx < logLimX(i)) cycle
            logPDF = logPDFNF(i) + (alpha(i) - 1._RKG) * logx
            return
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPiwiPowetoLogPDF_ENABLED && BAN_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0_IK < bin .and. bin <= size(alpha, kind = IK), SK_"@setPiwiPowetoLogPDFNF(): The conditions `0_IK < bin .and. bin <= size(alpha, kind = IK)` must hold. bin, size(alpha, kind = IK) = "//getStr([bin, size(alpha, kind = IK)]))
        CHECK_ASSERTION(__LINE__, size(alpha, kind = IK) == size(logPDFNF, kind = IK), SK_"@setPiwiPowetoLogPDFNF(): The condition `size(alpha) == size(logPDFNF)` must hold. size(alpha), size(logPDFNF) = "//getStr([size(alpha, kind = IK), size(logPDFNF, kind = IK)]))
        logPDF = logPDFNF(bin) + (alpha(bin) - 1._RKG) * logx

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPiwiPowetoCDF_ENABLED && ALDD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: cumSumArea(size(logLimX, kind = IK))
        call setPiwiPowetoCDF(cdf, logx, alpha, logLimX, getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea), cumSumArea)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getPiwiPowetoCDF_ENABLED && ALLC_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setPiwiPowetoCDF(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPiwiPowetoCDF_ENABLED && BAN_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_mathLogSubExp, only: getLogSubExp
        CHECK_ASSERTION(__LINE__, size(logLimX, kind = IK) == size(cumSumArea, kind = IK), SK_"@setPiwiPowetoCDF(): The condition `size(logLimX) == size(cumSumArea)` must hold. size(logLimX), size(cumSumArea) = "//getStr([size(logLimX, kind = IK), size(cumSumArea, kind = IK)]))
        CHECK_ASSERTION(__LINE__, size(logLimX, kind = IK) == size(alpha, kind = IK) + 1_IK, SK_"@setPiwiPowetoCDF(): The condition `size(logLimX) == size(alpha) + 1` must hold. size(logLimX), size(alpha) = "//getStr([size(logLimX, kind = IK), size(alpha, kind = IK)]))
        CHECK_ASSERTION(__LINE__, size(alpha, kind = IK) == size(logPDFNF, kind = IK), SK_"@setPiwiPowetoCDF(): The condition `size(alpha) == size(logPDFNF)` must hold. size(alpha), size(logPDFNF) = "//getStr([size(alpha, kind = IK), size(logPDFNF, kind = IK)]))
        CHECK_ASSERTION(__LINE__, 0_IK < bin .and. bin < size(logLimX, kind = IK), SK_"@setPiwiPowetoCDF(): The conditions `0_IK < bin .and. bin < size(logLimX)` must hold. bin, size(logLimX) = "//getStr([bin, size(logLimX, kind = IK)]))
        CHECK_ASSERTION(__LINE__, logLimX(bin) <= logx .and. logx <= logLimX(bin + 1), SK_"@setPiwiPowetoCDF(): The conditions `logLimX(bin) <= logx .and. logx < logLimX(bin + 1)` must hold. logLimX(bin), logx, logLimX(bin+1) = "//getStr([logLimX(bin), logx, logLimX(bin+1)]))
#if     CHECK_ENABLED
        block
            real(RKG) :: logPDFNF_ref(size(logPDFNF, kind = IK)), cumSumArea_ref(size(cumSumArea, kind = IK))
            logPDFNF_ref = getPiwiPowetoLogPDFNF(alpha, logLimX, cumSumArea_ref)
            CHECK_ASSERTION(__LINE__, all(abs(logPDFNF - logPDFNF_ref) <= epsilon(0._RKG) * 1000), SK_"@setPiwiPowetoCDF(): The condition `all(abs(logPDFNF - logPDFNF_ref) <= epsilon(0._RKG) * 1000)` must hold. logPDFNF - logPDFNF_ref = "//getStr([logPDFNF - logPDFNF_ref]))
            CHECK_ASSERTION(__LINE__, all(abs(cumSumArea - cumSumArea_ref) <= epsilon(0._RKG) * 1000), SK_"@setPiwiPowetoCDF(): The condition `all(abs(cumSumArea - cumSumArea_ref) <= epsilon(0._RKG) * 1000)` must hold. cumSumArea - cumSumArea_ref = "//getStr([cumSumArea - cumSumArea_ref]))
        end block
#endif
        cdf = cumSumArea(bin)
        if (logx /= logLimX(bin)) then
            if (alpha(bin) > 0._RKG) then
                cdf = cumSumArea(bin) + exp(logPDFNF(bin) + getLogSubExp(smaller = alpha(bin) * logLimX(bin), larger = alpha(bin) * logx)) / alpha(bin)
            elseif (alpha(bin) < 0._RKG) then
                cdf = cumSumArea(bin) - exp(logPDFNF(bin) + getLogSubExp(smaller = alpha(bin) * logx, larger = alpha(bin) * logLimX(bin))) / alpha(bin)
            else
                cdf = cumSumArea(bin) + exp(logPDFNF(bin) + log(logx - logLimX(bin)))
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setPiwiPowetoCDF_ENABLED && MAN_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenAlpha
        lenAlpha = size(alpha, kind = IK)
        CHECK_ASSERTION(__LINE__, lenAlpha > 0_IK, SK_"@setPiwiPowetoCDF(): The condition `size(alpha) > 0` must hold. size(logLimX) = "//getStr(lenAlpha))
        CHECK_ASSERTION(__LINE__, size(logPDFNF,1,IK) == lenAlpha, SK_"@setPiwiPowetoCDF(): The condition `size(logPDFNF,1,IK) == lenAlpha` must hold. size(logPDFNF,1,IK), lenAlpha = "//getStr([size(logPDFNF,1,IK), lenAlpha]))
        CHECK_ASSERTION(__LINE__, size(logLimX,1,IK) == lenAlpha + 1_IK, SK_"@setPiwiPowetoCDF(): The condition `size(logLimX) == size(alpha) + 1` must hold. size(logLimX,1,IK), size(alpha,1,IK) = "//getStr([size(logLimX,1,IK), size(alpha,1,IK)]))
        CHECK_ASSERTION(__LINE__, logLimX(1) <= logx .and. logx <= logLimX(lenAlpha + 1), SK_"@setPiwiPowetoCDF(): The condition `logLimX(1) <= logx and. logx < logLimX(size(logLimX))` must hold. logLimX(1), logx, logLimX(size(logLimX)) = "//getStr([logLimX(1), logx, logLimX(lenAlpha + 1)]))
        CHECK_ASSERTION(__LINE__, isAscending(logLimX), SK_"@setPiwiPowetoCDF(): The conditions `isAscending(logLimX)` must hold. logLimX = "//getStr(logLimX))
        !check_assertion(__LINE__, logLimX(size(logLimX,1,IK)) < log(huge(logLimX)), SK_"@getPiwiPowetoLogPDFNF(): The condition `logLimX(size(logLimX)) < log(huge(logLimX))` must hold. logLimX = "//getStr(logLimX))
        !check_assertion(__LINE__, merge(alpha(lenAlpha) < -1._RKG, .true., size(logLimX,1,IK) == lenAlpha), SK_"@getPiwiPowetoLogPDFNF(): The condition `alpha(size(alpha)) < -1._RKG` must hold. alpha = "//getStr(alpha))
        do i = lenAlpha, 1_IK, -1_IK
            if (logx < logLimX(i)) cycle
            call setPiwiPowetoCDF(cdf, logx, alpha, logLimX, logPDFNF, cumSumArea, i)
            return
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif