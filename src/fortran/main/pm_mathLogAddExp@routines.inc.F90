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
!>  \final
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
        complex(CKG), parameter :: ONE = cmplx(1._CKG, 0._CKG, CKG)
        complex(CKG)            :: ratio, y, z
#define GET_REAL(x) x%re
#elif   RK_ENABLED
        real(RKG)   , parameter :: ONE = 1._RKG
        real(RKG)               :: ratio, y, z
#define GET_REAL(x) x
#else
#error  "Unrecognized interface."
#endif
        integer     , parameter :: TKG = kind(logAddExp) ! type kind generic.
        real(TKG)   , parameter :: LOGTINY = log(tiny(0._TKG))
        real(TKG)   , parameter :: LOGEPS = log(epsilon(0._TKG))
        CHECK_ASSERTION(__LINE__, real(smaller, TKG) <= real(larger, TKG), \
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