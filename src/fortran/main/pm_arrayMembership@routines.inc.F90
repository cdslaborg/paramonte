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
!>  This file contains procedures implementations of the module [pm_arrayMembership](@ref pm_arrayMembership).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define size function.
#if     D0_D0_ENABLED
#define GET_SIZE(x) len(x, kind = IK)
#define GET_INDEX(i) i:i
#elif   D0_D1_ENABLED || D1_D1_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE(x) size(x, kind = IK)
#else
#error  "Unrecognized interface."
#endif
        ! Define equivalence check operator.
#if     LK_ENABLED
#define IS_NEQ .neqv.
#elif   SK_ENABLED || IK_ENABLED || CK_ENABLED || RK_ENABLED
#define IS_NEQ /=
#else
#error  "Unrecognized interface."
#endif
        ! Define comparison check operators.
#if     SK_ENABLED || IK_ENABLED || RK_ENABLED
#define IS_LESS(a, b) a < b
#define IS_MORE(a, b) a > b
#elif   LK_ENABLED
#define IS_LESS(a, b) .not. a .and. b
#define IS_MORE(a, b) a .and. .not. b
#elif   CK_ENABLED
#define IS_LESS(a, b) (a%re < b%re .or. (a%re == b%re .and. a%im < b%im))
#define IS_MORE(a, b) (a%re > b%re .or. (a%re == b%re .and. a%im > b%im))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     in_ENABLED && (D0_D0_ENABLED || D1_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, j, lenVal, lenSet
        lenVal = GET_SIZE(val)
        lenSet = GET_SIZE(set)
        loopVal: do i = 1, lenVal
            loopSet: do j = 1, lenSet
                if (set(GET_INDEX(j)) IS_NEQ val(GET_INDEX(i))) cycle loopSet
                member(i) = .true._LK
                cycle loopVal
            end do loopSet
            member(i) = .false._LK
        end do loopVal

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   in_ENABLED && D0_D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: j, lenSet
        lenSet = GET_SIZE(set)
        do j = 1, lenSet
            if (set(j) IS_NEQ val) cycle
            member = .true._LK
            return
        end do
        member = .false._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   inrange_ENABLED && (D0_D0_ENABLED || D1_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenVal
        lenVal = GET_SIZE(val)
        do i = 1, lenVal
            member(i) = .not. logical(IS_LESS(val(GET_INDEX(i)), set(GET_INDEX(1))) .or. IS_MORE(val(GET_INDEX(i)), set(GET_INDEX(2))), LK)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   inrange_ENABLED && D0_D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        member = .not. logical(IS_LESS(val, set(1)) .or. IS_MORE(val, set(2)), LK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   allin_ENABLED && (D0_D0_ENABLED || D1_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, j, lenVal, lenSet
        lenVal = GET_SIZE(val)
        lenSet = GET_SIZE(Set)
        loopVal: do i = 1, lenVal
            loopSet: do j = 1, lenSet
                if (Set(GET_INDEX(j)) IS_NEQ val(GET_INDEX(i))) cycle loopSet
                cycle loopVal
            end do loopSet
            allMember = .false._LK
            return
        end do loopVal
        allMember = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   allinrange_ENABLED && (D0_D0_ENABLED || D1_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenVal
        lenVal = GET_SIZE(val)
        do i = 1, lenVal
            if (IS_LESS(val(GET_INDEX(i)), Set(GET_INDEX(1))) .or. IS_MORE(val(GET_INDEX(i)), Set(GET_INDEX(2)))) then
                allMember = .false._LK
                return
            end if
        end do
        allMember = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   anyin_ENABLED && (D0_D0_ENABLED || D1_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, j, lenVal, lenSet
        lenVal = GET_SIZE(val)
        lenSet = GET_SIZE(Set)
        loopVal: do i = 1, lenVal
            loopSet: do j = 1, lenSet
                if (Set(GET_INDEX(j)) IS_NEQ val(GET_INDEX(i))) cycle loopSet
                anyMember = .true._LK
                return
            end do loopSet
        end do loopVal
        anyMember = .false._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   anyinrange_ENABLED && (D0_D0_ENABLED || D1_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: i, lenVal
        lenVal = GET_SIZE(val)
        do i = 1, lenVal
            if (IS_LESS(val(GET_INDEX(i)), Set(GET_INDEX(1))) .or. IS_MORE(val(GET_INDEX(i)), Set(GET_INDEX(2)))) cycle
            anyMember = .true._LK
            return
        end do
        anyMember = .false._LK

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_LESS
#undef  IS_MORE
#undef  IS_NEQ