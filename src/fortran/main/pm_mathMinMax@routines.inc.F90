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
!>  This include file contains procedure implementations of [pm_mathMinMax](@ref pm_mathMinMax).
!>
!>  \author
!>  \AmirShahmoradi, Sunday 3:33 AM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! DEfine the comparison operation.
#if     LK_ENABLED
#define IS_MORE(a,b) a .and. .not. b
#elif   CK_ENABLED
#define IS_MORE(a,b) a%re > b%re .or. (a%re == b%re .and. a%im > b%im)
#else
#define IS_MORE(a,b) a > b
#endif
        ! Define the pair.
#if     Pair_ENABLED
#define a pair(1)
#define b pair(2)
#elif   !Indi_ENABLED
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%
#if     getMinMax_ENABLED
        !%%%%%%%%%%%%%%%%

        if (IS_MORE(a, b)) then
            minMax(1) = b
            minMax(2) = a
        else
            minMax(1) = a
            minMax(2) = b
        end if

        !%%%%%%%%%%%%%%%%
#elif   setMinMax_ENABLED
        !%%%%%%%%%%%%%%%%

#if     SK_ENABLED
        character(max(len(a,IK),len(b,IK)),SKG) :: tmp
#elif   IK_ENABLED
        integer(IKG) :: tmp
#elif   LK_ENABLED
        logical(LKG) :: tmp
#elif   CK_ENABLED
        complex(CKG) :: tmp
#elif   RK_ENABLED
        real(RKG) :: tmp
#elif   !getMinMax_ENABLED
#error  "Unrecognized interface."
#endif
        if (IS_MORE(a, b)) then
            tmp = a
            a = b
            b = tmp
        end if
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  IS_MORE
#undef  a
#undef  b