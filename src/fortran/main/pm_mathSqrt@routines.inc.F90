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
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !getSqrt_ENABLED
#error  "Unrecognized interface."
#endif
#if     Bin_ENABLED || Def_ENABLED
        integer(IKG) :: middle, upper
        integer(IKG) , parameter :: maxIntSqrt = floor(sqrt(real(huge(0_IKG), RKB)), IKG)
        CHECK_ASSERTION(__LINE__, 0_IKG <= posint, SK_"@getSqrt(): The condition `0 <= posint` must hold. posint = "//getStr(posint))
        upper = min(maxIntSqrt, posint) + 1_IKG ! Avoid possible overflow.
        intSqrt = 0_IKG
        do
            if (upper - intSqrt == 1_IKG) return
            middle = (intSqrt + upper) / 2_IKG
            if (middle * middle <= posint) then
                intSqrt = middle
            else
                upper = middle
            end if
        end do
#elif   Lin_ENABLED
        ! linear search, ascending, using addition.
        ! (intSqrt + 1)^2 = intSqrt^2 + 2L + 1 = intSqrt^2 + 1 + \sum _{i = 1}^intSqrt 2
        integer(IKG) :: a, d
        integer(IKG) , parameter :: huge_IKG = huge(0_IKG)
        CHECK_ASSERTION(__LINE__, 0_IKG <= posint, SK_"@getSqrt(): The condition `0 <= posint` must hold. posint = "//getStr(posint))
        intSqrt = 0_IKG
        a = 1_IKG
        d = 3_IKG
        do
            if (posint < a) return
            if (huge_IKG - a < d) return
            a = a + d ! (a + 1_IKG) ** 2_IKG
            d = d + 2_IKG
            intSqrt = intSqrt + 1_IKG
        end do
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif