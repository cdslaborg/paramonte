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
!>  This include file contains procedure implementation of [pm_mathExp](@ref pm_mathExp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 25, 2015, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%
#if     isIntPow_ENABLED
        !%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0_IKG < absx, SK_"@getExpNext(): The condition `0 < absx` must hold. absx = "//getStr(absx))
#if     Def_ENABLED
        powisint = logical(popcnt(absx) == 1, LK)
#elif   Arb_ENABLED
        CHECK_ASSERTION(__LINE__, 1_IKG <= base, SK_"@getExpNext(): The condition `0 <= base` must hold. base = "//getStr(base))
        powisint = .false._LK
        block
            integer(IKG) :: quotient, dividend
            dividend = absx
            do
                quotient = dividend / base
                if (quotient * base /= dividend) return
                if (quotient == 1_IKG) exit
                dividend = quotient
            end do
            powisint = .true._LK
        end block
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getExpNext_ENABLED || getExpPrev_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the rounding mode.
#if     getExpNext_ENABLED
#define GET_ROUND(X) ceiling(X, kind = IKG)
#define EXPONENT expNext
#elif   getExpPrev_ENABLED
#define GET_ROUND(X) floor(X, kind = IKG)
#define EXPONENT expPrev
#else
#error  "Unrecognized interface."
#endif
        ! Compute the exponent.
#if     IK_ENABLED
        !>  \devnote
        !>  A `real` value of kind \RK32 can represent `integer` values as large as `huge(1_int128) = 170141183460469231731687303715884105727 = 1.70141183E+38 < huge(1._RK32) = 3.40282347E+38 << huge(1._RK64)`.<br>
        !>  One can envision a distant future human society with advanced computers capable of representing higher precision integer value for which \RK32 or \RK64 would be insufficient.<br>
        use pm_kind, only: RKG => RKH
        integer(IKG), parameter :: ZERO = 0_IKG
#define GET_REAL(x) real(x, RKG)
#elif   RK_ENABLED
#define GET_REAL(x) x
        integer, parameter :: IKG = IK
        real(RKG), parameter :: ZERO = 0._RKG
#else
#error  "Unrecognized interface."
#endif
        real(RKG), parameter :: INV_LOG_TWO = 1._RKG / log(2._RKG)
        CHECK_ASSERTION(__LINE__, 0 < absx, SK_"@getExpNext(): The condition `0 < absx` must hold. absx = "//getStr(absx))
        if (present(base)) then
            CHECK_ASSERTION(__LINE__, 1._RKG < GET_REAL(base), SK_"@getExpNext(): The condition `1._RKG < base` must hold. base = "//getStr(base))
            EXPONENT = GET_ROUND(log(GET_REAL(absx)) / log(GET_REAL(base)))
        else ! assume the base is 2
            EXPONENT = GET_ROUND(log(GET_REAL(absx)) * INV_LOG_TWO)
        end if
#undef  GET_ROUND
#undef  GET_REAL
#undef  EXPONENT

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif