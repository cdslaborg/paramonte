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
!>  This include file contains procedure implementations of [pm_mathLogSubExp](@ref pm_mathLogSubExp).
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
        use pm_math1mexp, only: get1mexp
#if     Sel_ENABLED && CK_ENABLED
        complex(CKG) :: logRatio
#define GET_REAL(x) x%re
#elif   Sel_ENABLED && RK_ENABLED
        real(RKG) :: logRatio
#define GET_REAL(x) x
#elif   !Seq_ENABLED
#error  "Unrecognized interface."
#endif
        integer, parameter :: TKG = kind(logSubExp) ! type kind generic.
        real(TKG), parameter :: LOGTINY = log(tiny(0._TKG))
        CHECK_ASSERTION(__LINE__, real(smaller, TKG) < real(larger, TKG), \
        SK_"@getLogSubExp(): The condition `real(smaller) < real(larger)` must hold. smaller, larger = "//getStr([smaller, larger]))
#if     Seq_ENABLED
        logSubExp = larger + log(get1mexp(smaller - larger)) ! log(1._TKG - exp(smaller - larger))
        !logSubExp = larger + log(1._TKG - exp(smaller - larger))
#elif   Sel_ENABLED
        logSubExp = larger
        logRatio = smaller - larger
        if (GET_REAL(logRatio) > LOGTINY) logSubExp = logSubExp + log(get1mexp(logRatio)) ! log(1._TKG - exp(logRatio))
#else
#error  "Unrecognized interface."
#endif
#undef  GET_REAL
#undef  smaller
#undef  larger