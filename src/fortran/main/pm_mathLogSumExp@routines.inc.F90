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
!>  This include file contains procedure implementation of
!>  [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp).
!>
!>  \author
!>  \AmirShahmoradi, Sunday 3:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define TYPE_KIND complex(TKC)
#elif   RK_ENABLED
#define TYPE_KIND real(TKC)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: i
        real(TKC), parameter :: LOGTINY = log(tiny(0._TKC))
#if     Sel_ENABLED
        TYPE_KIND :: logRatio
#elif   !Seq_ENABLED
#error  "Unrecognized interface."
#endif
#if     Def_ENABLED
        TYPE_KIND :: maxArray
        maxArray = maxval(real(array, TKC))
#elif   Max_ENABLED
        CHECK_ASSERTION(__LINE__, real(maxArray, TKC) == maxval(real(array, TKC)), SK_"@getLogSumExp(): The condition `real(maxArray) == maxval(real(array))` must hold. real(maxArray), maxval(real(array)) = "//getStr([real(maxArray, TKC), maxval(real(array, TKC))]))
#else
#error  "Unrecognized interface."
#endif
        logSumExp = 0._TKC
#if     Seq_ENABLED
        do i = 1_IK, size(array, kind = IK)
            logSumExp = logSumExp + exp(array(i) - maxArray)
        end do
        logSumExp = maxArray + log(logSumExp)
#elif   Sel_ENABLED
        do i = 1_IK, size(array, kind = IK)
            logRatio = array(i) - maxArray
            if (real(logRatio, kind = kind(logRatio)) > LOGTINY) logSumExp = logSumExp + exp(logRatio) ! fpp
        end do
        logSumExp = maxArray + log(logSumExp)
#else
#error  "Unrecognized interface."
#endif
#undef  TYPE_KIND