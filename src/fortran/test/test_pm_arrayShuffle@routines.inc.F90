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
!>  This file contains procedure implementations of [test_pm_arrayShuffle](@ref test_pm_arrayShuffle).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED
#define IS_NOT_EQUAL .neqv.
#else
#define IS_NOT_EQUAL /=
#endif

#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE len
#else
#define GET_SIZE size
#endif
        character(*, SK), parameter :: PROCEDURE_NAME = "@setShuffled()"
#if     SK_ENABLED && D0_ENABLED
#define ANY
        character(:,SKC), allocatable :: Array, arrayNew
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKC), dimension(:), allocatable :: Array, arrayNew
#elif   IK_ENABLED && D1_ENABLED
        integer(IKC)    , dimension(:), allocatable :: Array, arrayNew
#elif   LK_ENABLED && D1_ENABLED
        logical(LKC)    , dimension(:), allocatable :: Array, arrayNew
#elif   CK_ENABLED && D1_ENABLED
        complex(CKC)    , dimension(:), allocatable :: Array, arrayNew
#elif   RK_ENABLED && D1_ENABLED
        real(RKC)       , dimension(:), allocatable :: Array, arrayNew
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: count

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        call runTests(count)
        call runTests()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTests(count)
            integer(IK), intent(inout), optional :: count

            if (allocated(Array)) deallocate(Array) ! LCOV_EXCL_LINE

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            Array = ""
#elif       SK_ENABLED && D1_ENABLED
            allocate(character(2,SKC) :: Array(0))
#elif       IK_ENABLED && D1_ENABLED
            allocate(Array(0))
#elif       CK_ENABLED && D1_ENABLED
            allocate(Array(0))
#elif       RK_ENABLED && D1_ENABLED
            allocate(Array(0))
#elif       LK_ENABLED && D1_ENABLED
            allocate(Array(0))
#endif
            if (present(count)) count = 0_IK
            call report(count)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array has a shuffled array of length zero.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            Array = " "
#elif       SK_ENABLED && D1_ENABLED
            Array = [" "]
#elif       IK_ENABLED && D1_ENABLED
            Array = [1_IKC]
#elif       CK_ENABLED && D1_ENABLED
            Array = [(+1._CKC, -1._CKC)]
#elif       RK_ENABLED && D1_ENABLED
            Array = [1._RKC]
#elif       LK_ENABLED && D1_ENABLED
            Array = [.true._LKC]
#endif
            if (present(count)) call setUnifRand(count, 0_IK, GET_SIZE(Array, kind = IK))
            call report(count)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An array of length 1 has a shuffled array of length 1.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            Array = "ABCDEFGHIJK "
#elif       SK_ENABLED && D1_ENABLED
            Array = ["AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "II", "JJ", "KK", "  "]
#elif       IK_ENABLED && D1_ENABLED
            Array = [1_IKC, 2_IKC, 3_IKC, 4_IKC, 5_IKC, 6_IKC, 7_IKC, 8_IKC, 9_IKC]
#elif       CK_ENABLED && D1_ENABLED
            Array = [(+1._CKC, -1._CKC), (+2._CKC, -2._CKC), (+3._CKC, -3._CKC), (+4._CKC, -4._CKC), (+5._CKC, -5._CKC), (+6._CKC, -6._CKC), (+7._CKC, -7._CKC), (+8._CKC, -8._CKC), (+9._CKC, -9._CKC)]
#elif       RK_ENABLED && D1_ENABLED
            Array = [1._RKC, 2._RKC, 3._RKC, 4._RKC, 5._RKC, 6._RKC, 7._RKC, 8._RKC, 9._RKC]
#elif       LK_ENABLED && D1_ENABLED
            Array = [.false._LKC, .true._LKC, .false._LKC, .true._LKC, .false._LKC, .true._LKC, .false._LKC, .true._LKC, .false._LKC, .true._LKC, .false._LKC, .true._LKC]
#endif
            if (present(count)) call setUnifRand(count, 0_IK, GET_SIZE(Array, kind = IK))
            call report(count)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An array of arbitrary length must be shuffled properly.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(count)
            integer(IK), intent(in), optional :: count
            integer(IK) :: itry, count_def
#if         getShuffled_ENABLED
            itry = 0
            arrayNew = getShuffled(Array, count)
#elif       setShuffled_ENABLED
            type(xoshiro256ssw_type) :: rngx
            rngx = xoshiro256ssw_type()
            do itry = 1, 3
                arrayNew = Array
                if (itry == 1) then
                    call setShuffled(arrayNew, count)
                elseif  (itry == 2) then
                    call setShuffled(rngf_type(), arrayNew, count)
                elseif  (itry == 3) then
                    call setShuffled(rngx, arrayNew, count)
                end if
#else
#error          "Unrecognized interface."
#endif
                count_def = getOption(GET_SIZE(Array, kind = IK), count)
                assertion = assertion .and. (GET_SIZE(arrayNew) <= 1_IK .or. (arrayNew(1:count_def) .allin. Array))
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0,:,', '))")
                    write(test%disp%unit,"(*(g0,:,', '))") "itry       ", itry
                    write(test%disp%unit,"(*(g0,:,', '))") "Array      ", Array
                    write(test%disp%unit,"(*(g0,:,', '))") "arrayNew   ", arrayNew
                    write(test%disp%unit,"(*(g0,:,', '))")
                    ! LCOV_EXCL_STOP
                end if
#if         setShuffled_ENABLED
            end do
#endif
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  IS_NOT_EQUAL
#undef  GET_SIZE
#undef  ANY