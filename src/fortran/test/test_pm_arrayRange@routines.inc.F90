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
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED
#define ALL
#define TYPE_KIND character(1,SKG)
        integer(IK)                     :: step
        character(1,SKG)                :: start, finit
        character(:,SKG), allocatable   :: range, range_ref
#elif   IK_ENABLED
#define TYPE_KIND integer(IKG)
        integer(IKG)                    :: start, finit, step
        integer(IKG)    , allocatable   :: range(:), range_ref(:)
#elif   RK_ENABLED
#define TYPE_KIND real(RKG)
        real(RKG)                       :: start, finit, step
        real(RKG)       , allocatable   :: range(:), range_ref(:)
        real(RKG)       , parameter     :: TOL = 10 * epsilon(0._RKG)
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
        allocate(character(0,SKG) :: range_ref)
#elif   IK_ENABLED
        step = 1_IKG
        start = +0_IKG
        finit = -5_IKG
        allocate(range_ref(0))
#elif   RK_ENABLED
        step = 1._RKG
        start = +0._RKG
        finit = -5._RKG
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
        step = 1_IKG
        start = +1_IKG
        finit = +1_IKG
        range_ref = [start]
#elif   RK_ENABLED
        step = 1._RKG
        start = +1._RKG
        finit = +1._RKG
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
        step = 2_IKG
        start = +1_IKG
        finit = +2_IKG
        range_ref = [start]
#elif   RK_ENABLED
        step = 2._RKG
        start = +1._RKG
        finit = +2._RKG
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
        range_ref = SKG_"aceg"
#elif   IK_ENABLED
        step = 2_IKG
        start = -1_IKG
        finit = +5_IKG
        range_ref = [-1_IKG, +1_IKG, +3_IKG, +5_IKG]
#elif   RK_ENABLED
        step = 2._RKG
        start = -1._RKG
        finit = +6._RKG
        range_ref = [real(RKG) :: -1, +1, +3, +5]
#endif
        call report(start, finit, step = step)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit, step = "//getStr([start, finit])//SK_", "//getStr(step))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        step = -3
        start = "h"
        finit = "a"
        range_ref = SKG_"heb"
#elif   IK_ENABLED
        step = -3_IKG
        start = +6_IKG
        finit = -1_IKG
        range_ref = [6_IKG, +3_IKG, 0_IKG]
#elif   RK_ENABLED
        step = -3._RKG
        start = +6._RKG
        finit = -1._RKG
        range_ref = [real(RKG) :: 6, +3, 0]
#endif
        call report(start, finit, step = step)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit, step = "//getStr([start, finit])//SK_", "//getStr(step))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        start = "c"
        finit = "a"
        range_ref = SKG_"cba"
#elif   IK_ENABLED
        start = +3_IKG
        finit = -1_IKG
        range_ref = [integer(IKG) :: 3, 2, 1, 0, -1]
#elif   RK_ENABLED
        start = 1._RKG
        finit = nearest(nearest(1._RKG, -1._RKG), -1._RKG)
        range_ref = [real(RKG) :: start, nearest(1._RKG, -1._RKG), finit]
#endif
        call report(start, finit)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit = "//getStr([start, finit]))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        start = "A"
        finit = "z"
#elif   IK_ENABLED
        start = 1_IKG
        finit = 10_IKG  
#elif   RK_ENABLED
        start = 1._RKG
        finit = 1._RKG + 1000 * epsilon(0._RKG)
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
        range_ref = SKG_"abc"
#elif   IK_ENABLED
        step = 1_IKG
        start = +1_IKG
        finit = +3_IKG
        range_ref = [integer(IKG) :: 1, 2, 3]
#elif   RK_ENABLED
        step = 1._RKG
        start = 1._RKG
        finit = 4._RKG
        range_ref = [real(RKG) :: 1, 2, 3]
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
        range_ref = SKG_"fdb"
#elif   IK_ENABLED
        step = -2_IKG
        start = +6_IKG
        finit = +1_IKG
        range_ref = [integer(IKG) :: 6, 4, 2]
#elif   RK_ENABLED
        step = -2._RKG
        start = 6._RKG
        finit = 1._RKG
        range_ref = [real(RKG) :: 6, 4, 2]
#endif
        allocate(range, mold = range_ref)
        call report(start, step)
        call test%assert(assertion, SK_"getRange() must yield a proper `range` with start, finit, step = "//getStr([start, finit])//SK_", "//getStr(step))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED
        start = "a"
        range_ref = SKG_"abc"
#elif   IK_ENABLED
        start = +1_IKG
        range_ref = [integer(IKG) :: 1, 2, 3]
#elif   RK_ENABLED
        start = 1._RKG
        allocate(range_ref(3))
        block
            integer(IK) :: i
            range_ref(1) = nearest(1._RKG, 1._RKG)
            do i = 2, size(range_ref)
                range_ref(i) = nearest(range_ref(i - 1), 1._RKG)
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