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
!>  This module contains implementations of the tests of the procedures under the generic interface [setRemapped](@ref pm_arrayRemap::setRemapped).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define comparison operator.
#if     LK_ENABLED && D1_ENABLED
#define ISEQ .eqv.
#else
#define ISEQ ==
#endif
        ! Define sizing function.
#if     !(SK_ENABLED && D0_ENABLED)
#define GET_SIZE size
#endif
        ! Define procedure name.
#if     getRemapped_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@getRemapped()"
#elif   setRemapped_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@setRemapped()"
#endif
#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE len
#define ALL
        character(:,SKC), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKC), dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   IK_ENABLED && D1_ENABLED
        integer(IKC)    , dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   LK_ENABLED && D1_ENABLED
        logical(LKC)    , dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   CK_ENABLED && D1_ENABLED
        complex(CKC)    , dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   RK_ENABLED && D1_ENABLED
        real(RKC)       , dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#else
#error  "Unrecognized interface."
#endif
        integer(IK), allocatable    :: index(:)
        logical(LK)                 :: backward_def
        logical(LK)                 :: arrayNewEnabled

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        arrayNewEnabled = .false._LK
        call runTestsWith()
        call runTestsWith(action = reverse)
        arrayNewEnabled = .true._LK
        call runTestsWith()
        call runTestsWith(action = reverse)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(action)
            type(reverse_type), intent(in), optional :: action

            backward_def = present(action)
            if (allocated(index)) deallocate(index)
            if (allocated(Array)) deallocate(Array)
            if (allocated(ArrayNew_ref)) deallocate(ArrayNew_ref)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            Array = ""
            ArrayNew_ref = ""
#elif       SK_ENABLED && D1_ENABLED
            allocate(character(2,SKC) :: Array(0), ArrayNew_ref(0))
#elif       IK_ENABLED && D1_ENABLED
            allocate(Array(0), ArrayNew_ref(0))
#elif       LK_ENABLED && D1_ENABLED
            allocate(Array(0), ArrayNew_ref(0))
#elif       CK_ENABLED && D1_ENABLED
            allocate(Array(0), ArrayNew_ref(0))
#elif       RK_ENABLED && D1_ENABLED
            allocate(Array(0), ArrayNew_ref(0))
#endif
            allocate(index(GET_SIZE(Array)))
            call report()
            call test%assert(assertion, desc = PROCEDURE_NAME//SK_": An empty array has a remapped array of length zero.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            Array = SKC_" "
            ArrayNew_ref = SKC_" "
#elif       SK_ENABLED && D1_ENABLED
            Array = [SKC_" "]
            ArrayNew_ref = [SKC_" "]
#elif       IK_ENABLED && D1_ENABLED
            Array = [1_IKC]
            ArrayNew_ref = [1_IKC]
#elif       LK_ENABLED && D1_ENABLED
            Array = [.true._LKC]
            ArrayNew_ref = [.true._LKC]
#elif       CK_ENABLED && D1_ENABLED
            Array = [(+1._CKC, -1._CKC)]
            ArrayNew_ref = [(+1._CKC, -1._CKC)]
#elif       RK_ENABLED && D1_ENABLED
            Array = [1._RKC]
            ArrayNew_ref = [1._RKC]
#endif
            index = [1_IK]
            call report(action)
            call test%assert(assertion, desc = PROCEDURE_NAME//SK_": An array of length 1 has a remapped array of length 1.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            Array = SKC_"ABCDE "
            if (backward_def) then
                ArrayNew_ref = SKC_" EDCBA"
            else
                ArrayNew_ref = Array
            end if
#elif       SK_ENABLED && D1_ENABLED
            Array = ["AA", "BB", "CC", "DD", "EE", "  "]
            if (backward_def) then
                ArrayNew_ref = Array(size(Array):1:-1)
            else
                ArrayNew_ref = Array
            end if
#elif       IK_ENABLED && D1_ENABLED
            Array = [1_IKC, 2_IKC, 3_IKC, 4_IKC, 5_IKC, 6_IKC]
            if (backward_def) then
                ArrayNew_ref = Array(size(Array):1:-1)
            else
                ArrayNew_ref = Array
            end if
#elif       LK_ENABLED && D1_ENABLED
            Array = [.false._LKC, .true._LKC, .false._LKC, .true._LKC, .false._LKC, .true._LKC]
            if (backward_def) then
                ArrayNew_ref = Array(size(Array):1:-1)
            else
                ArrayNew_ref = Array
            end if
#elif       CK_ENABLED && D1_ENABLED
            Array = [(+1._CKC, -1._CKC), (+2._CKC, -2._CKC), (+3._CKC, -3._CKC), (+4._CKC, -4._CKC), (+5._CKC, -5._CKC), (+6._CKC, -6._CKC)]
            if (backward_def) then
                ArrayNew_ref = Array(size(Array):1:-1)
            else
                ArrayNew_ref = Array
            end if
#elif       RK_ENABLED && D1_ENABLED
            Array = [1._RKC, 2._RKC, 3._RKC, 4._RKC, 5._RKC, 6._RKC]
            if (backward_def) then
                ArrayNew_ref = Array(size(Array):1:-1)
            else
                ArrayNew_ref = Array
            end if
#endif
            index = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK]
            call report(action)
            call test%assert(assertion, desc = PROCEDURE_NAME//SK_": The order of an array of length 6 must not change by a map that is the same as the array indices, unless `action = reverse`.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            Array = "ABCDE "
            if (backward_def) then
                ArrayNew_ref = "ADEC B"
            else
                ArrayNew_ref = "B CEDA"
            end if
#elif       SK_ENABLED && D1_ENABLED
            Array = ["AA", "BB", "CC", "DD", "EE", "  "]
            if (backward_def) then
                ArrayNew_ref = ["AA", "DD", "EE", "CC", "  ", "BB"]
            else
                ArrayNew_ref = ["BB", "  ", "CC", "EE", "DD", "AA"]
            end if
#elif       IK_ENABLED && D1_ENABLED
            Array = [1_IKC, 2_IKC, 3_IKC, 4_IKC, 5_IKC, 6_IKC]
            if (backward_def) then
                ArrayNew_ref = [1_IKC, 4_IKC, 5_IKC, 3_IKC, 6_IKC, 2_IKC]
            else
                ArrayNew_ref = [2_IKC, 6_IKC, 3_IKC, 5_IKC, 4_IKC, 1_IKC]
            end if
#elif       LK_ENABLED && D1_ENABLED
            Array = [.false._LKC, .true._LKC, .false._LKC, .true._LKC, .false._LKC, .true._LKC]
            if (backward_def) then
                ArrayNew_ref = [.false._LKC, .true._LKC, .false._LKC, .false._LKC, .true._LKC, .true._LKC]
            else
                ArrayNew_ref = [.true._LKC, .true._LKC, .false._LKC, .false._LKC, .true._LKC, .false._LKC]
            end if
#elif       CK_ENABLED && D1_ENABLED
            Array = [(+1._CKC, -1._CKC), (+2._CKC, -2._CKC), (+3._CKC, -3._CKC), (+4._CKC, -4._CKC), (+5._CKC, -5._CKC), (+6._CKC, -6._CKC)]
            if (backward_def) then
                ArrayNew_ref = [(+1._CKC, -1._CKC), (+4._CKC, -4._CKC), (+5._CKC, -5._CKC), (+3._CKC, -3._CKC), (+6._CKC, -6._CKC), (+2._CKC, -2._CKC)]
            else
                ArrayNew_ref = [(+2._CKC, -2._CKC), (+6._CKC, -6._CKC), (+3._CKC, -3._CKC), (+5._CKC, -5._CKC), (+4._CKC, -4._CKC), (+1._CKC, -1._CKC)]
            end if
#elif       RK_ENABLED && D1_ENABLED
            Array = [1._RKC, 2._RKC, 3._RKC, 4._RKC, 5._RKC, 6._RKC]
            if (backward_def) then
                ArrayNew_ref = [1._RKC, 4._RKC, 5._RKC, 3._RKC, 6._RKC, 2._RKC]
            else
                ArrayNew_ref = [2._RKC, 6._RKC, 3._RKC, 5._RKC, 4._RKC, 1._RKC]
            end if
#endif
            index = [2_IK, 6_IK, 3_IK, 5_IK, 4_IK, 1_IK]
            call report(action)
            call test%assert(assertion, desc = PROCEDURE_NAME//SK_": An array of length 6 must be remapped correctly with unique indices.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            Array = "ABCDE "
            if (backward_def) then
                ArrayNew_ref = "AAAAAC"
            else
                ArrayNew_ref = "CAAAAA"
            end if
#elif       SK_ENABLED && D1_ENABLED
            Array = ["AA", "BB", "CC", "DD", "EE", "  "]
            if (backward_def) then
                ArrayNew_ref = ["AA", "AA", "AA", "AA", "AA", "CC"]
            else
                ArrayNew_ref = ["CC", "AA", "AA", "AA", "AA", "AA"]
            end if
#elif       IK_ENABLED && D1_ENABLED
            Array = [1_IKC, 2_IKC, 3_IKC, 4_IKC, 5_IKC, 6_IKC]
            if (backward_def) then
                ArrayNew_ref = [1_IKC, 1_IKC, 1_IKC, 1_IKC, 1_IKC, 3_IKC]
            else
                ArrayNew_ref = [3_IKC, 1_IKC, 1_IKC, 1_IKC, 1_IKC, 1_IKC]
            end if
#elif       LK_ENABLED && D1_ENABLED
            Array = [.false._LKC, .true._LKC, .false._LKC, .true._LKC, .false._LKC, .true._LKC]
            if (backward_def) then
                ArrayNew_ref = [.false._LKC, .false._LKC, .false._LKC, .false._LKC, .false._LKC, .false._LKC]
            else
                ArrayNew_ref = [.false._LKC, .false._LKC, .false._LKC, .false._LKC, .false._LKC, .false._LKC]
            end if
#elif       CK_ENABLED && D1_ENABLED
            Array = [(+1._CKC, -1._CKC), (+2._CKC, -2._CKC), (+3._CKC, -3._CKC), (+4._CKC, -4._CKC), (+5._CKC, -5._CKC), (+6._CKC, -6._CKC)]
            if (backward_def) then
                ArrayNew_ref = [(+1._CKC, -1._CKC), (+1._CKC, -1._CKC), (+1._CKC, -1._CKC), (+1._CKC, -1._CKC), (+1._CKC, -1._CKC), (+3._CKC, -3._CKC)]
            else
                ArrayNew_ref = [(+3._CKC, -3._CKC), (+1._CKC, -1._CKC), (+1._CKC, -1._CKC), (+1._CKC, -1._CKC), (+1._CKC, -1._CKC), (+1._CKC, -1._CKC)]
            end if
#elif       RK_ENABLED && D1_ENABLED
            Array = [1._RKC, 2._RKC, 3._RKC, 4._RKC, 5._RKC, 6._RKC]
            if (backward_def) then
                ArrayNew_ref = [1._RKC, 1._RKC, 1._RKC, 1._RKC, 1._RKC, 3._RKC]
            else
                ArrayNew_ref = [3._RKC, 1._RKC, 1._RKC, 1._RKC, 1._RKC, 1._RKC]
            end if
#endif
            index = [3_IK, 1_IK, 1_IK, 1_IK, 1_IK, 1_IK]
            call report(action)
            call test%assert(assertion, desc = PROCEDURE_NAME//SK_": An array of length 6 must be remapped correctly with unique indices.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(action)

            type(reverse_type), intent(in), optional :: action
            if (allocated(arrayNew)) deallocate(arrayNew)
            allocate(arrayNew, source = Array)

            if (present(action) .and. arrayNewEnabled) then
#if             getRemapped_ENABLED
                arrayNew = getRemapped(Array, index, action)
#elif           setRemapped_ENABLED
                call setRemapped(Array, index, action, arrayNew)
#else
#error          "Unrecognized interface."
#endif
            elseif (present(action)) then
#if             getRemapped_ENABLED
                arrayNew = getRemapped(arrayNew, index, action)
#elif           setRemapped_ENABLED
                call setRemapped(arrayNew, index, action)
#endif
            elseif (arrayNewEnabled) then
#if             getRemapped_ENABLED
                arrayNew = getRemapped(Array, index)
#elif           setRemapped_ENABLED
                call setRemapped(Array, index, arrayNew = arrayNew)
#endif
            else
#if             getRemapped_ENABLED
                arrayNew = getRemapped(arrayNew, index)
#elif           setRemapped_ENABLED
                call setRemapped(arrayNew, index)
#endif
            end if

            assertion = assertion .and. ALL(arrayNew ISEQ ArrayNew_ref)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Array              ", Array
                write(test%disp%unit,"(*(g0,:,', '))") "arrayNew           ", arrayNew
                write(test%disp%unit,"(*(g0,:,', '))") "ArrayNew_ref       ", ArrayNew_ref
                write(test%disp%unit,"(*(g0,:,', '))") "backward_def       ", backward_def
                write(test%disp%unit,"(*(g0,:,', '))") "arrayNewEnabled    ", arrayNewEnabled
                write(test%disp%unit,"(*(g0,:,', '))") "index              ", index
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  GET_SIZE
#undef  ISEQ
#undef  ALL
