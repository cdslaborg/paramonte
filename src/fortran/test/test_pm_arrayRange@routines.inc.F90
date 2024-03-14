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
!>  This include file contains the implementations of the tests of procedures with generic interfaces [pm_arrayRange](@ref pm_arrayRange).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED
#define ALL
#define TYPE_KIND character(1,SKC)
        integer(IK)                     :: step
        character(1,SKC)                :: start, finit
        character(:,SKC), allocatable   :: range, range_ref
#elif   IK_ENABLED
#define TYPE_KIND integer(IKC)
        integer(IKC)                    :: start, finit, step
        integer(IKC)    , allocatable   :: range(:), range_ref(:)
#elif   RK_ENABLED
#define TYPE_KIND real(RKC)
        real(RKC)                       :: start, finit, step
        real(RKC)       , allocatable   :: range(:), range_ref(:)
        real(RKC)       , parameter     :: TOL = 10 * epsilon(0._RKC)
#else
#error  "Unrecognized interface."
#endif
        assertion = .true._LK

        !%%%%%%%%%%%%%%%
#if     getRange_ENABLED
        !%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        step = 1
        start = "z"
        finit = "a"
        allocate(character(0,SKC) :: range_ref)
#elif   IK_ENABLED
        step = 1_IKC
        start = +0_IKC
        finit = -5_IKC
        allocate(range_ref(0))
#elif   RK_ENABLED
        step = 1._RKC
        start = +0._RKC
        finit = -5._RKC
        allocate(range_ref(0))
#endif
        call report(start, finit, step = step)
        call test%assert(assertion, SK_"getRange() must yield an empty `range` with `finit < start` with a positive `step`.")

        call report(finit, start, step = -step)
        call test%assert(assertion, SK_"getRange() must yield an empty `range` with `finit > start` with a negative `step`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        step = 1
        start = "a"
        finit = "a"
        range_ref = start
#elif   IK_ENABLED
        step = 1_IKC
        start = +1_IKC
        finit = +1_IKC
        range_ref = [start]
#elif   RK_ENABLED
        step = 1._RKC
        start = +1._RKC
        finit = +1._RKC
        range_ref = [start]
#endif

        call report(start, finit)
        call test%assert(assertion, SK_"getRange() must yield an `range` of size `1` with `finit == start`.")

        call report(start, finit, step = step)
        call test%assert(assertion, SK_"getRange() must yield an `range` of size `1` with `finit == start` with an equality `step`.")

        call report(finit, start, step = -step)
        call test%assert(assertion, SK_"getRange() must yield an `range` of size `1` with `finit == start` with an equality `step`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        step = 2
        start = "a"
        finit = "b"
        range_ref = start
#elif   IK_ENABLED
        step = 2_IKC
        start = +1_IKC
        finit = +2_IKC
        range_ref = [start]
#elif   RK_ENABLED
        step = 2._RKC
        start = +1._RKC
        finit = +2._RKC
        range_ref = [start]
#endif

        call report(start, finit, step = step)
        call test%assert(assertion, SK_"getRange() must yield an `range` of size `1` with `finit < start + step`.")

        range_ref = finit
        call report(finit, start, step = -step)
        call test%assert(assertion, SK_"getRange() must yield an `range` of size `1` with `finit - step < start`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        step = 2
        start = "a"
        finit = "h"
        range_ref = SKC_"aceg"
#elif   IK_ENABLED
        step = 2_IKC
        start = -1_IKC
        finit = +5_IKC
        range_ref = [-1_IKC, +1_IKC, +3_IKC, +5_IKC]
#elif   RK_ENABLED
        step = 2._RKC
        start = -1._RKC
        finit = +6._RKC
        range_ref = [real(RKC) :: -1, +1, +3, +5]
#endif
        call report(start, finit, step = step)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit, step = "//getStr([start, finit])//SK_", "//getStr(step))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        step = -3
        start = "h"
        finit = "a"
        range_ref = SKC_"heb"
#elif   IK_ENABLED
        step = -3_IKC
        start = +6_IKC
        finit = -1_IKC
        range_ref = [6_IKC, +3_IKC, 0_IKC]
#elif   RK_ENABLED
        step = -3._RKC
        start = +6._RKC
        finit = -1._RKC
        range_ref = [real(RKC) :: 6, +3, 0]
#endif
        call report(start, finit, step = step)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit, step = "//getStr([start, finit])//SK_", "//getStr(step))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        start = "c"
        finit = "a"
        range_ref = SKC_"cba"
#elif   IK_ENABLED
        start = +3_IKC
        finit = -1_IKC
        range_ref = [integer(IKC) :: 3, 2, 1, 0, -1]
#elif   RK_ENABLED
        start = 1._RKC
        finit = nearest(nearest(1._RKC, -1._RKC), -1._RKC)
        range_ref = [real(RKC) :: start, nearest(1._RKC, -1._RKC), finit]
#endif
        call report(start, finit)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit = "//getStr([start, finit]))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        start = "A"
        finit = "z"
#elif   IK_ENABLED
        start = 1_IKC
        finit = 10_IKC  
#elif   RK_ENABLED
        start = 1._RKC
        finit = 1._RKC + 1000 * epsilon(0._RKC)
#endif
        range = getRange(start, finit)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit = "//getStr([start, finit]))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(start, finit, step)
            TYPE_KIND, intent(in) :: start, finit
#if         SK_ENABLED
            integer(IK), intent(in), optional :: step
#else
            TYPE_KIND, intent(in), optional  :: step
#endif
            type(display_type) :: disp
            if (present(step)) then
                range = getRange(start, finit, step)
            else
                range = getRange(start, finit)
            end if
#if         SK_ENABLED || IK_ENABLED
            assertion = assertion .and. ALL(range == range_ref)
#elif       RK_ENABLED
            assertion = assertion .and. all(isClose(range_ref, range, abstol = TOL))
#endif
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip
                call disp%show("start")
                call disp%show( start )
                call disp%show("finit")
                call disp%show( finit )
                call disp%show("present(step)")
                call disp%show( present(step) )
                if (present(step)) then
                call disp%show("step")
                call disp%show( step )
                end if
                call disp%show("range")
                call disp%show( range )
                call disp%show("range_ref")
                call disp%show( range_ref )
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%
#elif   setRange_ENABLED
        !%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        step = 1_IK
        start = "a"
        finit = "c"
        range_ref = SKC_"abc"
#elif   IK_ENABLED
        step = 1_IKC
        start = +1_IKC
        finit = +3_IKC
        range_ref = [integer(IKC) :: 1, 2, 3]
#elif   RK_ENABLED
        step = 1._RKC
        start = 1._RKC
        finit = 4._RKC
        range_ref = [real(RKC) :: 1, 2, 3]
#endif
        allocate(range, mold = range_ref)
        call report(start, step)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit, step = "//getStr([start, finit])//SK_", "//getStr(step))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        step = -2_IK
        start = "f"
        finit = "a"
        range_ref = SKC_"fdb"
#elif   IK_ENABLED
        step = -2_IKC
        start = +6_IKC
        finit = +1_IKC
        range_ref = [integer(IKC) :: 6, 4, 2]
#elif   RK_ENABLED
        step = -2._RKC
        start = 6._RKC
        finit = 1._RKC
        range_ref = [real(RKC) :: 6, 4, 2]
#endif
        allocate(range, mold = range_ref)
        call report(start, step)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit, step = "//getStr([start, finit])//SK_", "//getStr(step))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        start = "a"
        range_ref = SKC_"abc"
#elif   IK_ENABLED
        start = +1_IKC
        range_ref = [integer(IKC) :: 1, 2, 3]
#elif   RK_ENABLED
        start = 1._RKC
        allocate(range_ref(3))
        block
            integer(IK) :: i
            range_ref(1) = nearest(1._RKC, 1._RKC)
            do i = 2, size(range_ref)
                range_ref(i) = nearest(range_ref(i - 1), 1._RKC)
            end do
        end block
#endif
        allocate(range, mold = range_ref)
        call report(start)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit, step = "//getStr([start, finit])//SK_", "//getStr(step))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(start, step)
            TYPE_KIND, intent(in) :: start
#if         SK_ENABLED
            integer(IK), intent(in), optional :: step
#else
            TYPE_KIND, intent(in), optional  :: step
#endif
            type(display_type) :: disp
            if (present(step)) then
                call setRange(range, start, step)
            else
                call setRange(range, start)
            end if
#if         SK_ENABLED || IK_ENABLED
            assertion = assertion .and. ALL(range == range_ref)
#elif       RK_ENABLED
            assertion = assertion .and. all(isClose(range_ref, range, abstol = TOL))
#endif
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip
                call disp%show("start")
                call disp%show( start )
                call disp%show("finit")
                call disp%show( finit )
                call disp%show("present(step)")
                call disp%show( present(step) )
                if (present(step)) then
                call disp%show("step")
                call disp%show( step )
                end if
                call disp%show("range")
                call disp%show( range )
                call disp%show("range_ref")
                call disp%show( range_ref )
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else
#error  "Unrecognized interface."
#endif

        subroutine reset()
            if (allocated(range)) deallocate(range)
            if (allocated(range_ref)) deallocate(range_ref)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  TYPE_KIND
#undef  ALL