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
!>  This include file contains implementations of the procedures in module [pm_mathSqrt](@ref pm_mathSqrt).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !getSqrt_ENABLED
#error  "Unrecognized interface."
#endif
#if     Bin_ENABLED || Def_ENABLED
        integer(IKC) :: middle, upper
        integer(IKC) , parameter :: maxIntSqrt = floor(sqrt(real(huge(0_IKC), RKB)), IKC)
        CHECK_ASSERTION(__LINE__, 0_IKC <= posint, SK_"@getSqrt(): The condition `0 <= posint` must hold. posint = "//getStr(posint))
        upper = min(maxIntSqrt, posint) + 1_IKC ! Avoid possible overflow.
        intSqrt = 0_IKC
        do
            if (upper - intSqrt == 1_IKC) return
            middle = (intSqrt + upper) / 2_IKC
            if (middle * middle <= posint) then
                intSqrt = middle
            else
                upper = middle
            end if
        end do
#elif   Lin_ENABLED
        ! linear search, ascending, using addition.
        ! (intSqrt + 1)^2 = intSqrt^2 + 2L + 1 = intSqrt^2 + 1 + \sum _{i = 1}^intSqrt 2
        integer(IKC) :: a, d
        integer(IKC) , parameter :: huge_IKC = huge(0_IKC)
        CHECK_ASSERTION(__LINE__, 0_IKC <= posint, SK_"@getSqrt(): The condition `0 <= posint` must hold. posint = "//getStr(posint))
        intSqrt = 0_IKC
        a = 1_IKC
        d = 3_IKC
        do
            if (posint < a) return
            if (huge_IKC - a < d) return
            a = a + d ! (a + 1_IKC) ** 2_IKC
            d = d + 2_IKC
            intSqrt = intSqrt + 1_IKC
        end do
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif