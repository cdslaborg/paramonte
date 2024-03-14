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
!>  This include file contains procedure implementations of [pm_mathLogAddExp](@ref pm_mathLogAddExp).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     MM_ENABLED
#define smaller minMax(1)
#define larger minMax(2)
#elif   !SL_ENABLED
#error  "Unrecognized interface."
#endif
#if     CK_ENABLED
        complex(CKC), parameter :: ONE = cmplx(1._CKC, 0._CKC, CKC)
        complex(CKC)            :: ratio, y, z
#define GET_REAL(x) x%re
#elif   RK_ENABLED
        real(RKC)   , parameter :: ONE = 1._RKC
        real(RKC)               :: ratio, y, z
#define GET_REAL(x) x
#else
#error  "Unrecognized interface."
#endif
        integer     , parameter :: TKC = kind(logAddExp) ! This kind current.
        real(TKC)   , parameter :: LOGTINY = log(tiny(0._TKC))
        real(TKC)   , parameter :: LOGEPS = log(epsilon(0._TKC))
        CHECK_ASSERTION(__LINE__, real(smaller, TKC) <= real(larger, TKC), \
        SK_"@getLogAddExp(): The condition `real(smaller) <= real(larger)` must hold. smaller, larger = "//getStr([smaller, larger]))
#if     Seq_ENABLED
        ratio = exp(smaller - larger)
        y = ONE + ratio
        z = y - ONE
        logAddExp = larger + log(y) - (z - ratio) / y
#elif   Sel_ENABLED
        logAddExp = larger
        ratio = smaller - larger
        if (GET_REAL(ratio) > LOGEPS) then
            logAddExp = logAddExp + log(ONE + exp(ratio))
        elseif (GET_REAL(ratio) > LOGTINY) then
            ratio = exp(ratio)
            y = ONE + ratio
            z = y - ONE
            logAddExp = logAddExp + log(y) - (z - ratio) / y
        end if
#else
#error  "Unrecognized interface."
#endif

#undef  GET_REAL
#undef  smaller
#undef  larger