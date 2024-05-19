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
!>  This include file contains procedure implementations of the tests of [pm_mathCumSum](@ref pm_mathCumSum).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: itry
#if     getCumSum_ENABLED
        logical(LK), parameter :: isnew = .true.
#elif   setCumSum_ENABLED
        logical(LK) :: isnew
#else
#error  "Unrecognized interface."
#endif
#if     IK_ENABLED
#define TYPE_KIND integer(TKG)
        integer(TKG), parameter :: TOL = 0_TKG
        integer(TKG), parameter :: LB = -1_TKG, UB = 1_TKG, ZERO = 0_TKG
        integer(TKG), allocatable :: cumSum_ref(:), cumSum(:), array(:), diff(:)
#elif   CK_ENABLED
#define TYPE_KIND complex(TKG)
        complex(TKG), parameter :: TOL = epsilon(1._TKG) * 10._TKG * (1._TKG, 1._TKG)
        complex(TKG), parameter :: LB = (-1._TKG, -2._TKG), UB = (2._TKG, 1._TKG), ZERO = (0._TKG, 0._TKG)
        complex(TKG), allocatable :: cumSum_ref(:), cumSum(:), array(:), diff(:)
#elif   RK_ENABLED
#define TYPE_KIND real(TKG)
        real(TKG), parameter :: TOL = epsilon(1._TKG) * 10._TKG
        real(TKG), parameter :: LB = -1._TKG, UB = 1._TKG, ZERO = 0._TKG
        real(TKG), allocatable :: cumSum_ref(:), cumSum(:), array(:), diff(:)
#else
#error  "Unrecognized interface."
#endif
        logical(LK) :: isbackward, isreverse

        assertion = .true._LK

        do itry = 1, 300

#if         setCumSum_ENABLED
            isnew = getUnifRand()
#endif
            isreverse = getUnifRand()
            isbackward = getUnifRand()
            array = getUnifRand(LB, UB, getUnifRand(1_IK, 10_IK))
            call setResized(cumSum, size(array, 1, IK))

            if (isbackward .and. isreverse) then
                cumSum_ref = getCumSum_ref(array, backward, reverse)
#if             setCumSum_ENABLED
                if (isnew) then
                    call setCumSum(cumSum, array, backward, reverse)
                else
                    cumSum = array
                    call setCumSum(cumSum, array, backward, reverse)
                end if
#elif           getCumSum_ENABLED
                cumSum = getCumSum(array, backward, reverse)
#else
#error          "Unrecognized interface."
#endif
            elseif (isbackward) then
                cumSum_ref = getCumSum_ref(array, backward, nothing)
#if             setCumSum_ENABLED
                if (isnew) then
                    call setCumSum(cumSum, array, backward, nothing)
                else
                    cumSum = array
                    call setCumSum(cumSum, backward, nothing)
                end if
#elif           getCumSum_ENABLED
                cumSum = getCumSum(array, backward, nothing)
#endif
            elseif (isreverse) then
                cumSum_ref = getCumSum_ref(array, forward, reverse)
#if             setCumSum_ENABLED
                if (isnew) then
                    call setCumSum(cumSum, array, forward, reverse)
                else
                    cumSum = array
                    call setCumSum(cumSum, forward, reverse)
                end if
#elif           getCumSum_ENABLED
                cumSum = getCumSum(array, forward, reverse)
#endif
            else
                cumSum_ref = getCumSum_ref(array, forward, nothing)
#if             setCumSum_ENABLED
                if (isnew) then
                    call setCumSum(cumSum, array, forward, nothing)
                else
                    cumSum = array
                    call setCumSum(cumSum, forward, nothing)
                end if
#elif           getCumSum_ENABLED
                cumSum = getCumSum(array, forward, nothing)
#endif
            end if
            call report(__LINE__)

#if         getCumSum_ENABLED
            if (isbackward .and. .not. isreverse) call runTestsWith(direction = backward)
            if (isreverse .and. .not. isbackward) call runTestsWith(action = reverse)
            if (.not. (isbackward .or. isreverse)) call runTestsWith()
#endif

        end do

    contains

#if     getCumSum_ENABLED
        subroutine runTestsWith(direction, action)
            class(action_type), intent(in), optional :: action
            class(direction_type), intent(in), optional :: direction
            cumSum = getCumSum(array, direction, action)
            call report(__LINE__)
        end subroutine
#endif

        function getCumSum_ref(array, direction, action) result(cumSum)
            class(direction_type), intent(in), optional :: direction
            class(action_type), intent(in), optional :: action
            class(direction_type), allocatable :: direction_def
            class(action_type), allocatable :: action_def
            TYPE_KIND, intent(in) :: array(:)
            TYPE_KIND :: cumSum(size(array, 1, IK))
            integer(IK) :: i
            action_def = nothing
            direction_def = forward
            if (present(action)) action_def = action
            if (present(direction)) direction_def = direction
            cumSum = array
            if (same_type_as(direction_def, backward)) call setReversed(cumSum)
            do i = 2, size(array, 1, IK)
                cumSum(i) = cumSum(i) + cumSum(i - 1)
            end do
            if (same_type_as(action_def, reverse)) call setReversed(cumSum)
        end function

        subroutine report(line)
            integer :: line
            diff = cumSum - cumSum_ref
            assertion = all(-TOL <= diff .and. diff <= TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip
                call test%disp%show("isnew")
                call test%disp%show( isnew )
                call test%disp%show("isreverse")
                call test%disp%show( isreverse )
                call test%disp%show("isbackward")
                call test%disp%show( isbackward )
                call test%disp%show("cumSum_ref")
                call test%disp%show( cumSum_ref )
                call test%disp%show("cumSum")
                call test%disp%show( cumSum )
                call test%disp%show("array")
                call test%disp%show( array )
                call test%disp%show("diff")
                call test%disp%show( diff )
                call test%disp%show("TOL")
                call test%disp%show( TOL )
                call test%disp%skip
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The output `cumSum` must be correctly computed.", line)
        end subroutine
#undef  TYPE_KIND