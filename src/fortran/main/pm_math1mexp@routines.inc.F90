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
!>  This include file contains procedure implementations of [pm_math1mexp](@ref pm_math1mexp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
        complex(CKG), parameter :: ONE = cmplx(1._CKG, 0._CKG, CKG), ZERO = cmplx(0._CKG, 0._CKG, CKG)
#define GET_REAL(x) x%re
#elif   RK_ENABLED
        real(RKG)   , parameter :: ONE = 1._RKG, ZERO = 0._RKG
#define GET_REAL(x) x
#else
#error  "Unrecognized interface."
#endif
        integer , parameter :: TKG = kind(onemexp) ! type kind generic.
#if     Seq_ENABLED && CK_ENABLED
        complex(CKG) :: tsterm, i
#elif   Seq_ENABLED && RK_ENABLED
        real(RKG) :: tsterm, i
#elif   Sel_ENABLED
        real(TKG), parameter :: NEG_LOG_HUGE = -log(huge(0._TKG))
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, real(x, TKG) < log(huge(0._TKG)), \
        SK_"@get1mexp(): The condition `real(x, TKG) <= huge(0._TKG)` must hold. x = "//getStr(x))
#if     Seq_ENABLED
        if (abs(GET_REAL(x)) < log(2._TKG)) then
            onemexp = x
            tsterm = x
            i = 1._TKG
            do
                i = i + ONE
                tsterm = tsterm * x / i
                onemexp = onemexp + tsterm
                if (abs(GET_REAL(tsterm)) > abs(GET_REAL(onemexp)) * epsilon(0._TKG)) cycle
                exit
            end do
            onemexp = -onemexp
        else
            onemexp = ONE - exp(x)
        end if
#elif   Sel_ENABLED
        ! Is this really needed? any number smaller than tiny? Yes: zero
        if (abs(GET_REAL(x)) < tiny(0._TKG)) then
            onemexp = ONE
        elseif (GET_REAL(x) < NEG_LOG_HUGE) then
            onemexp = ZERO
        else
            onemexp = get1mexp(x)
        end if
#else
#error  "Unrecognized interface."
#endif

#undef  GET_REAL