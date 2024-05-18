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
!>  This include file contains procedure implementations of the tests of [pm_mathCumPropExp](@ref pm_mathCumPropExp).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: itry
#if     getCumPropExp_ENABLED
        logical(LK), parameter :: isnew = .true.
#elif   setCumPropExp_ENABLED
        logical(LK) :: isnew
#else
#error  "Unrecognized interface."
#endif
#if     Sel_ENABLED
        type(selection_type), parameter :: control = selection_type()
#elif   Seq_ENABLED
        type(sequence_type), parameter :: control = sequence_type()
#else
#error  "Unrecognized interface."
#endif
        real(RKC), parameter :: TOL = epsilon(1._RKC) * 10._RKC
        !real(RKC), parameter :: LB = 1._RKC, UB = 2._RKC
        real(RKC), parameter :: LB = minexponent(0._RKC), UB = maxexponent(0._RKC)
        real(RKC), allocatable :: cumPropExp_ref(:), cumPropExp(:), array(:), diff(:)
        logical(LK) :: isbackward, isreverse

        assertion = .true._LK

        do itry = 1, 300

#if         setCumPropExp_ENABLED
            isnew = getUnifRand()
#endif
            isreverse = getUnifRand()
            isbackward = getUnifRand()
            array = getUnifRand(LB, UB, getUnifRand(1_IK, 10_IK))
            call setResized(cumPropExp, size(array, 1, IK))

            if (isbackward .and. isreverse) then
                cumPropExp_ref = getCumPropExp_ref(array, maxval(array), backward, reverse)
#if             setCumPropExp_ENABLED
                if (isnew) then
                    call setCumPropExp(cumPropExp, array, maxval(array), control, backward, reverse)
                else
                    cumPropExp = array
                    call setCumPropExp(cumPropExp, array, maxval(array), control, backward, reverse)
                end if
#elif           getCumPropExp_ENABLED
                cumPropExp = getCumPropExp(array, maxval(array), control, backward, reverse)
#else
#error          "Unrecognized interface."
#endif
            elseif (isbackward) then
                cumPropExp_ref = getCumPropExp_ref(array, maxval(array), backward, nothing)
#if             setCumPropExp_ENABLED
                if (isnew) then
                    call setCumPropExp(cumPropExp, array, maxval(array), control, backward, nothing)
                else
                    cumPropExp = array
                    call setCumPropExp(cumPropExp, maxval(array), control, backward, nothing)
                end if
#elif           getCumPropExp_ENABLED
                cumPropExp = getCumPropExp(array, maxval(array), control, backward, nothing)
#endif
            elseif (isreverse) then
                cumPropExp_ref = getCumPropExp_ref(array, maxval(array), forward, reverse)
#if             setCumPropExp_ENABLED
                if (isnew) then
                    call setCumPropExp(cumPropExp, array, maxval(array), control, forward, reverse)
                else
                    cumPropExp = array
                    call setCumPropExp(cumPropExp, maxval(array), control, forward, reverse)
                end if
#elif           getCumPropExp_ENABLED
                cumPropExp = getCumPropExp(array, maxval(array), control, forward, reverse)
#endif
            else
                cumPropExp_ref = getCumPropExp_ref(array, maxval(array), forward, nothing)
#if             setCumPropExp_ENABLED
                if (isnew) then
                    call setCumPropExp(cumPropExp, array, maxval(array), control, forward, nothing)
                else
                    cumPropExp = array
                    call setCumPropExp(cumPropExp, maxval(array), control, forward, nothing)
                end if
#elif           getCumPropExp_ENABLED
                cumPropExp = getCumPropExp(array, maxval(array), control, forward, nothing)
#endif
            end if
            call report(__LINE__)

#if         getCumPropExp_ENABLED
            if (isbackward .and. .not. isreverse) call runTestsWith(direction = backward)
            if (isreverse .and. .not. isbackward) call runTestsWith(action = reverse)
            if (.not. (isbackward .or. isreverse)) call runTestsWith()
#endif

        end do

    contains

#if     getCumPropExp_ENABLED
        subroutine runTestsWith(direction, action)
            class(action_type), intent(in), optional :: action
            class(direction_type), intent(in), optional :: direction
            cumPropExp = getCumPropExp(array, maxval(array), control, direction, action)
            call report(__LINE__)
            cumPropExp = getCumPropExp(array, maxval(array), direction, action)
            call report(__LINE__)
        end subroutine
#endif

        function getCumPropExp_ref(array, maxArray, direction, action) result(cumPropExp)
            class(direction_type), intent(in), optional :: direction
            class(action_type), intent(in), optional :: action
            class(direction_type), allocatable :: direction_def
            class(action_type), allocatable :: action_def
            real(RKC), intent(in) :: array(:), maxArray
            real(RKC) :: cumPropExp(size(array, 1, IK))
            action_def = nothing
            direction_def = forward
            if (present(action)) action_def = action
            if (present(direction)) direction_def = direction
            if (same_type_as(direction_def, backward) .and. same_type_as(action_def, reverse)) then
                cumPropExp = getCumSum(exp(array - maxArray), backward, reverse)
            elseif (same_type_as(direction_def, backward) .and. same_type_as(action_def, nothing)) then
                cumPropExp = getCumSum(exp(array - maxArray), backward, nothing)
            elseif (same_type_as(direction_def, forward) .and. same_type_as(action_def, reverse)) then
                cumPropExp = getCumSum(exp(array - maxArray), forward, reverse)
            elseif (same_type_as(direction_def, forward) .and. same_type_as(action_def, nothing)) then
                cumPropExp = getCumSum(exp(array - maxArray), forward, nothing)
            end if
            cumPropExp = cumPropExp / maxval(cumPropExp)
        end function

        subroutine report(line)
            integer :: line
            diff = abs(cumPropExp - cumPropExp_ref)
            assertion = all(diff < TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip
                call test%disp%show("isnew")
                call test%disp%show( isnew )
                call test%disp%show("isreverse")
                call test%disp%show( isreverse )
                call test%disp%show("isbackward")
                call test%disp%show( isbackward )
                call test%disp%show("cumPropExp_ref")
                call test%disp%show( cumPropExp_ref )
                call test%disp%show("cumPropExp")
                call test%disp%show( cumPropExp )
                call test%disp%show("array")
                call test%disp%show( array )
                call test%disp%show("diff")
                call test%disp%show( diff )
                call test%disp%show("TOL")
                call test%disp%show( TOL )
                call test%disp%skip
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The output `cumPropExp` must be correctly computed.", line)
        end subroutine