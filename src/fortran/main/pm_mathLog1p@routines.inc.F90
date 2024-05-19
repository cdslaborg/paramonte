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
!>  This include file contains procedure implementations of [pm_mathLog1p](@ref pm_mathLog1p).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     Seq_ENABLED && CK_ENABLED
        complex(CKG) :: y, z
#elif   Seq_ENABLED && RK_ENABLED
        real(RKG) :: y, z
!#elif   Sel_ENABLED && CK_ENABLED
!        use pm_complexCompareAll, only: operator(<)
!        complex(CKG) :: absX
!#elif   Sel_ENABLED && RK_ENABLED
!        real(RKG) :: absX
#else
#error  "Unrecognized interface."
#endif
#if     CK_ENABLED
        complex(CKG), parameter :: ONE = cmplx(1._CKG, 0._CKG, CKG)
        complex(CKG), parameter :: ZERO = cmplx(0._CKG, 0._CKG, CKG)
        complex(CKG), parameter :: EPSX = cmplx(epsilon(0._CKG), epsilon(0._CKG), CKG)
        complex(CKG), parameter :: TINYX = cmplx(tiny(0._CKG), tiny(0._CKG), CKG)
!#define GET_REAL(x) x%re
#elif   RK_ENABLED
        real(RKG)   , parameter :: ONE = 1._RKG, ZERO = 0._RKG, EPSX = epsilon(x), TINYX = tiny(x)
!#define GET_REAL(x) x
#else
#error  "Unrecognized interface."
#endif
        integer , parameter :: TKG = kind(x) ! type kind generic.
        CHECK_ASSERTION(__LINE__, real(x, TKG) > -1._TKG, SK_"@getLog1p(): The condition `real(x) > -1.` must hold. x = "//getStr(x))
        CHECK_ASSERTION(__LINE__, real(x, TKG) < huge(0._TKG), SK_"@getLog1p(): The condition `real(x) <= huge(x)` must hold. x = "//getStr(x))
#if     Seq_ENABLED
        y = ONE + x
        z = y - ONE
        log1p = log(y) - (z - x) / y
!#elif   Sel_ENABLED
!        !absRealX = abs(GET_REAL(x))
!        absX = abs(x)
!        ! Is this really needed? any number smaller than tiny? Yes: zero
!        ! if (absRealX < tiny(0._TKG)) then
!        if (absX < TINYX) then
!            log1p = ZERO
!        !elseif (absRealX < epsilon(0._TKG)) then
!        elseif (absX < EPSX) then
!            log1p = getLog1p(x)
!        else
!            log1p = log(ONE + x)
!        end if
#else
#error  "Unrecognized interface."
#endif

!#undef  GET_REAL