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
!>  This file contains procedure implementations of [pm_arrayRange](@ref pm_arrayRange).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 PM, September 22, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     !(getRange_ENABLED || setRange_ENABLED)
#error  "Unrecognized interface."
#endif
        ! Set the sizing function.
#if     D0_ENABLED && SK_ENABLED
        use pm_kind, only: IKC => IK
        integer(IKC) :: index
#define GET_INDEX(i) i:i
#define GET_SIZE len
#elif   D1_ENABLED && (IK_ENABLED || RK_ENABLED)
#define GET_INDEX(i) i
#define GET_SIZE size
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     (SK_ENABLED || IK_ENABLED) && (D0_ENABLED || D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IKC) :: lenRange, irange
#if     Unit_ENABLED && getRange_ENABLED
        integer(IKC) :: step
        if (start < stop) then
            step = 1_IKC
        else
            step = -1_IKC
        end if
#elif   Unit_ENABLED
        integer(IKC), parameter :: step = 1_IKC
#elif   Step_ENABLED
        CHECK_ASSERTION(__LINE__, step /= 0_IKC, SK_"@setRange(): The condition `step /= 0` must hold. step = "//getStr(step))
#else
#error  "Unrecognized interface."
#endif
        lenRange = GET_SIZE(range, kind = IKC)
#if     SK_ENABLED && D0_ENABLED
        index = ichar(start, IKC)
#endif
        if (0_IKC < lenRange) then
            range(GET_INDEX(1_IKC)) = start
            do irange = 2_IKC, lenRange
#if             SK_ENABLED
                index = index + step
                range(GET_INDEX(irange)) = char(index, SKC)
#elif           IK_ENABLED
                range(GET_INDEX(irange)) = range(GET_INDEX(irange - 1_IKC)) + step
#else
#error          "Unrecognized interface."
#endif
            end do
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   D1_ENABLED && RK_ENABLED && Unit_ENABLED && getRange_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKC) :: direction, next
        integer(IK) :: iell, nell
        direction = stop - start
        if (0._RKC < direction) then
            iell = 1; nell = 127
            call setResized(range, nell)
            range(iell) = start
            do
                if (iell < nell) then
                    next = nearest(range(iell), direction)
                    if (stop < next) exit
                    iell = iell + 1
                    range(iell) = next
                else
                    nell = 2 * nell
                    call setResized(range, nell)
                end if
            end do
            range = range(1:iell)
        elseif (direction < 0._RKC) then
            iell = 1; nell = 127
            call setResized(range, nell)
            range(iell) = start
            do
                if (iell < nell) then
                    next = nearest(range(iell), direction)
                    if (next < stop) exit
                    iell = iell + 1
                    range(iell) = next
                else
                    nell = 2 * nell
                    call setResized(range, nell)
                end if
            end do
            range = range(1:iell)
        elseif (0._RKC == direction) then
            range = [start]
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   D1_ENABLED && RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        real(RKC), parameter :: direction = 1._RKC
        integer(IK) :: iell, nell
#if     Step_ENABLED
        CHECK_ASSERTION(__LINE__, step /= 0._RKC, SK_"@setRange(): The condition `step /= 0` must hold. step = "//getStr(step))
#endif
        nell = size(range, 1, IK)
        if (nell == 0_IK) return
        range(1) = start
        do iell = 2, nell
#if         Unit_ENABLED && setRange_ENABLED
            range(iell) = nearest(range(iell - 1), direction)
#elif       Step_ENABLED && (getRange_ENABLED || setRange_ENABLED)
            range(iell) = range(iell - 1) + step
#else
#error      "Unrecognized interface."
#endif
        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_INDEX
#undef  GET_SIZE