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
!>  This include file contains implementations of the procedures in module [pm_mathUnsigned](@ref pm_mathUnsigned).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     uadd_ENABLED && IK_ENABLED
        integer(IKG) :: aa, bb
        integer(IKG), parameter :: HUGE_IKG = huge(0_IKG) !, SIGNBITLOC = bit_size(a) - 1_IKG
        CHECK_ASSERTION(__LINE__, 0_IKG <= a, SK_"@getFactoring(): The condition `0 <= a` must hold for corresponding input arguments. a = "//getStr(a))
        CHECK_ASSERTION(__LINE__, 0_IKG <= b, SK_"@getFactoring(): The condition `0 <= b` must hold for corresponding input arguments. b = "//getStr(b))
        if (HUGE_IKG - a < b) then ! overflowed.
            !aa = ibset(a, SIGNBITLOC)
            !bb = ibset(b, SIGNBITLOC)
            aa = a - HUGE_IKG - 1_IKG
            bb = b - HUGE_IKG - 1_IKG
            sum = aa + bb
        else
            sum = a + b
        end if
#else
#error  "Unrecognized interface."
#endif