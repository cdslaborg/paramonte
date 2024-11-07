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
!>  This module contains implementations of the tests of the procedures under the generic interfaces<br>
!>  [setReversed](@ref pm_arrayReverse::setReversed)<br>
!>  [getReversed](@ref pm_arrayReverse::getReversed)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getReversed_D1_LK_ENABLED || setReversed_D1_LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     !(getReversed_D0_SK_ENABLED || setReversed_D0_SK_ENABLED)
#define GET_SIZE size
#endif

#if     getReversed_D0_SK_ENABLED || setReversed_D0_SK_ENABLED
#define GET_SIZE len
#define ALL
        character(:,SKG), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   getReversed_D1_SK_ENABLED || setReversed_D1_SK_ENABLED
        character(2,SKG), dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   getReversed_D1_IK_ENABLED || setReversed_D1_IK_ENABLED
        integer(IKG)    , dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   getReversed_D1_CK_ENABLED || setReversed_D1_CK_ENABLED
        complex(CKG)    , dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   getReversed_D1_RK_ENABLED || setReversed_D1_RK_ENABLED
        real(RKG)       , dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#elif   getReversed_D1_LK_ENABLED || setReversed_D1_LK_ENABLED
        logical(LKG)    , dimension(:), allocatable :: Array, arrayNew, ArrayNew_ref
#else
#error  "Unrecognized interface."
#endif

#if     getReversed_D1_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@getReversed()"
#elif   setReversed_D1_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@setReversed()"
        logical(LK) :: arrayNewEnabled
        arrayNewEnabled = .false._LK
#else
#error  "Unrecognized Interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        call runTestsWith()
#if     setReversed_D1_ENABLED
        arrayNewEnabled = .true._LK
        call runTestsWith()
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith()

            if (allocated(Array)) deallocate(Array)
            if (allocated(ArrayNew_ref)) deallocate(ArrayNew_ref)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getReversed_D0_SK_ENABLED || setReversed_D0_SK_ENABLED
            Array = ""
            ArrayNew_ref = ""
#elif       getReversed_D1_SK_ENABLED || setReversed_D1_SK_ENABLED
            allocate(character(2,SKG) :: Array(0), ArrayNew_ref(0))
#elif       getReversed_D1_IK_ENABLED || setReversed_D1_IK_ENABLED
            allocate(Array(0), ArrayNew_ref(0))
#elif       getReversed_D1_CK_ENABLED || setReversed_D1_CK_ENABLED
            allocate(Array(0), ArrayNew_ref(0))
#elif       getReversed_D1_RK_ENABLED || setReversed_D1_RK_ENABLED
            allocate(Array(0), ArrayNew_ref(0))
#elif       getReversed_D1_LK_ENABLED || setReversed_D1_LK_ENABLED
            allocate(Array(0), ArrayNew_ref(0))
#endif
            call report()
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array has a reversed array of length zero.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getReversed_D0_SK_ENABLED || setReversed_D0_SK_ENABLED
            Array = " "
            ArrayNew_ref = " "
#elif       getReversed_D1_SK_ENABLED || setReversed_D1_SK_ENABLED
            Array = [" "]
            ArrayNew_ref = [" "]
#elif       getReversed_D1_IK_ENABLED || setReversed_D1_IK_ENABLED
            Array = [1_IKG]
            ArrayNew_ref = [1_IKG]
#elif       getReversed_D1_CK_ENABLED || setReversed_D1_CK_ENABLED
            Array = [(+1._CKG, -1._CKG)]
            ArrayNew_ref = [(+1._CKG, -1._CKG)]
#elif       getReversed_D1_RK_ENABLED || setReversed_D1_RK_ENABLED
            Array = [1._RKG]
            ArrayNew_ref = [1._RKG]
#elif       getReversed_D1_LK_ENABLED || setReversed_D1_LK_ENABLED
            Array = [.true._LKG]
            ArrayNew_ref = [.true._LKG]
#endif
            call report()
            call test%assert(assertion, PROCEDURE_NAME//SK_": An array of length 1 has a reversed array of length 1.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getReversed_D0_SK_ENABLED || setReversed_D0_SK_ENABLED
            Array = "ABCDE "
            ArrayNew_ref = " EDCBA"
#elif       getReversed_D1_SK_ENABLED || setReversed_D1_SK_ENABLED
            Array = ["AA", "BB", "CC", "DD", "EE", "  "]
            ArrayNew_ref = Array(size(Array):1:-1)
#elif       getReversed_D1_IK_ENABLED || setReversed_D1_IK_ENABLED
            Array = [1_IKG, 2_IKG, 3_IKG, 4_IKG, 5_IKG, 6_IKG]
            ArrayNew_ref = Array(size(Array):1:-1)
#elif       getReversed_D1_CK_ENABLED || setReversed_D1_CK_ENABLED
            Array = [(+1._CKG, -1._CKG), (+2._CKG, -2._CKG), (+3._CKG, -3._CKG), (+4._CKG, -4._CKG), (+5._CKG, -5._CKG), (+6._CKG, -6._CKG)]
            ArrayNew_ref = Array(size(Array):1:-1)
#elif       getReversed_D1_RK_ENABLED || setReversed_D1_RK_ENABLED
            Array = [1._RKG, 2._RKG, 3._RKG, 4._RKG, 5._RKG, 6._RKG]
            ArrayNew_ref = Array(size(Array):1:-1)
#elif       getReversed_D1_LK_ENABLED || setReversed_D1_LK_ENABLED
            Array = [.false._LKG, .true._LKG, .false._LKG, .true._LKG, .false._LKG, .true._LKG]
            ArrayNew_ref = Array(size(Array):1:-1)
#endif
            call report()
            call test%assert(assertion, PROCEDURE_NAME//SK_": An ordered array must be correctly reversed.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()

#if         getReversed_D1_ENABLED
            arrayNew = getReversed(Array)
#elif       setReversed_D1_ENABLED
            if (arrayNewEnabled) then
                if (allocated(arrayNew)) deallocate(arrayNew)
                allocate(arrayNew, mold = Array)
                call setReversed(Array, arrayNew)
            else
                if (allocated(arrayNew)) deallocate(arrayNew)
                allocate(arrayNew, source = Array)
                call setReversed(arrayNew)
            end if
#else
#error      "Unrecognized Interface."
#endif

            assertion = assertion .and. ALL(arrayNew IS_EQUAL ArrayNew_ref)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Array              ", Array
                write(test%disp%unit,"(*(g0,:,', '))") "arrayNew           ", arrayNew
                write(test%disp%unit,"(*(g0,:,', '))") "ArrayNew_ref       ", ArrayNew_ref
#if             setReversed_D1_ENABLED
                write(test%disp%unit,"(*(g0,:,', '))") "arrayNewEnabled    ", arrayNewEnabled
#endif
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  IS_EQUAL
#undef  GET_SIZE
#undef  ALL
