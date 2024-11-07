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
!>  This module contains implementations of the tests of the procedures under the generic interface
!>  [setInserted](@ref pm_arrayInsert::setInserted).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setInserted_D0_SK_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#else
#define GET_INDEX(i) i
#define GET_SIZE size
#endif

#if     setInserted_D1_LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

        use pm_val2str, only: getStr
        use pm_kind, only: LK, SK

        character(*, SK), parameter                 :: PROCEDURE_NAME = "@setInserted()"

#if     setInserted_D0_SK_ENABLED
#define ALL
        character(:,SKG), allocatable               :: array, insertion, arrayNew, arrayNewS_ref, arrayNewV_ref
#elif   setInserted_D1_SK_ENABLED
        character(2,SKG), allocatable, dimension(:) :: array, insertion, arrayNew, arrayNewS_ref, arrayNewV_ref
#elif   setInserted_D1_IK_ENABLED
        integer(IKG)    , allocatable, dimension(:) :: array, insertion, arrayNew, arrayNewS_ref, arrayNewV_ref
#elif   setInserted_D1_CK_ENABLED
        complex(CKG)    , allocatable, dimension(:) :: array, insertion, arrayNew, arrayNewS_ref, arrayNewV_ref
#elif   setInserted_D1_RK_ENABLED
        real(RKG)       , allocatable, dimension(:) :: array, insertion, arrayNew, arrayNewS_ref, arrayNewV_ref
#elif   setInserted_D1_LK_ENABLED
        logical(LKG)    , allocatable, dimension(:) :: array, insertion, arrayNew, arrayNewS_ref, arrayNewV_ref
#else
#error  "Unrecognized interface."
#endif
        integer(IK)     , allocatable               :: index(:)
        logical(LK)                                 :: getInsertedEnabled
        integer(IK) :: i

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK

        getInsertedEnabled = .false._LK
        do i = 1, 2
            call runInsertionTestsWith()
            call runInsertionTestsWith(sorted = .true._LK)
            call runInsertionTestsWith(sorted = .false._LK)
            call runInsertionTestsWith(positive = .true._LK)
            call runInsertionTestsWith(positive = .false._LK)
            call runInsertionTestsWith(positive = .true._LK, sorted = .true._LK)
            call runInsertionTestsWith(positive = .false._LK, sorted = .true._LK)
            call runInsertionTestsWith(positive = .false._LK, sorted = .false._LK)
            call runInsertionTestsWith(positive = .true._LK, sorted = .false._LK)
            getInsertedEnabled = .not. getInsertedEnabled
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runInsertionTestsWith(positive, sorted)

            use pm_option, only: getOption
            logical(LK), intent(in), optional :: positive, sorted

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         setInserted_D0_SK_ENABLED
            insertion = " "
            allocate(character(0,SKG) :: arrayNewV_ref, array)
#elif       setInserted_D1_SK_ENABLED
            insertion = [" "]
            allocate(character(2,SKG) :: arrayNewV_ref(0), array(0))
#elif       setInserted_D1_IK_ENABLED
            insertion = [1_IKG]
            allocate(arrayNewV_ref(0), array(0))
#elif       setInserted_D1_CK_ENABLED
            insertion = [1._CKG]
            allocate(arrayNewV_ref(0), array(0))
#elif       setInserted_D1_RK_ENABLED
            insertion = [1._RKG]
            allocate(arrayNewV_ref(0), array(0))
#elif       setInserted_D1_LK_ENABLED
            insertion = [.false._LKG]
            allocate(arrayNewV_ref(0), array(0))
#endif
            allocate(index(0))
            arrayNewS_ref = arrayNewV_ref

            call runTestWith(positive, sorted)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `arrayNew` with vector `insertion` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            call runTestWith(positive, sorted, scalarInsertionIndex = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `arrayNew` with scalar `insertion` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         setInserted_D0_SK_ENABLED
            allocate(character(0,SKG) :: arrayNewV_ref, array, insertion)
#elif       setInserted_D1_SK_ENABLED
            allocate(character(2,SKG) :: arrayNewV_ref(0), array(0), insertion(0))
#elif       setInserted_D1_IK_ENABLED
            allocate(arrayNewV_ref(0), array(0), insertion(0))
#elif       setInserted_D1_CK_ENABLED
            allocate(arrayNewV_ref(0), array(0), insertion(0))
#elif       setInserted_D1_RK_ENABLED
            allocate(arrayNewV_ref(0), array(0), insertion(0))
#elif       setInserted_D1_LK_ENABLED
            allocate(arrayNewV_ref(0), array(0), insertion(0))
#endif
            allocate(index(0))
            arrayNewS_ref = arrayNewV_ref

            call runTestWith(positive, sorted)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `arrayNew` with vector `insertion` of length zero with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         setInserted_D0_SK_ENABLED
            array = "AAAA"
            insertion = "X"
#elif       setInserted_D1_SK_ENABLED
            array = ["AA", "AA"]
            insertion = ["XX"]
#elif       setInserted_D1_IK_ENABLED
            array = [1_IKG, 1_IKG]
            insertion = [2_IKG]
#elif       setInserted_D1_CK_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            insertion = [(2._CKG,-2._CKG)]
#elif       setInserted_D1_RK_ENABLED
            array = [1._RKG, 1._RKG]
            insertion = [2._RKG]
#elif       setInserted_D1_LK_ENABLED
            array = [.false._LK, .false._LK]
            insertion = [.true._LK]
#endif
            allocate(index(0))
            arrayNewS_ref = array
            arrayNewV_ref = array

            call runTestWith(positive, sorted)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a single-element vector `insertion` with vector `insertion` with an empty `index` must yield an `arrayNew` that is identical to `array` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            call runTestWith(positive, sorted, scalarInsertionIndex = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a single-element scalar `insertion` with vector `insertion` with an empty `index` must yield an `arrayNew` that is identical to `array` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         setInserted_D0_SK_ENABLED
            array = "AA"
            insertion = "XY"
            arrayNewS_ref = "XAXA"
            arrayNewV_ref = "XYAXYA"
#elif       setInserted_D1_SK_ENABLED
            array = ["AA", "AA"]
            insertion = ["XX", "YY"]
            arrayNewS_ref = ["XX", "AA", "XX", "AA"]
            arrayNewV_ref = ["XX", "YY", "AA", "XX", "YY", "AA"]
#elif       setInserted_D1_IK_ENABLED
            array = [1_IKG, 1_IKG]
            insertion = [2_IKG, 3_IKG]
            arrayNewS_ref = [2_IKG, 1_IKG, 2_IKG, 1_IKG]
            arrayNewV_ref = [2_IKG, 3_IKG, 1_IKG, 2_IKG, 3_IKG, 1_IKG]
#elif       setInserted_D1_CK_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            insertion = [(2._CKG,-2._CKG), (3._CKG,-3._CKG)]
            arrayNewS_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG)]
            arrayNewV_ref = [(2._CKG,-2._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG)]
#elif       setInserted_D1_RK_ENABLED
            array = [1._RKG, 1._RKG]
            insertion = [2._RKG, 3._RKG]
            arrayNewS_ref = [2._RKG, 1._RKG, 2._RKG, 1._RKG]
            arrayNewV_ref = [2._RKG, 3._RKG, 1._RKG, 2._RKG, 3._RKG, 1._RKG]
#elif       setInserted_D1_LK_ENABLED
            array = [.false._LK, .false._LK]
            insertion = [.true._LK, .true._LK]
            arrayNewS_ref = [.true._LK, .false._LK, .true._LK, .false._LK]
            arrayNewV_ref = [.true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK]
#endif
            index = [1_IK, 2_IK]

            call runTestWith(positive, sorted)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a double-element vector `insertion` with `index = [1,2]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            call runTestWith(positive, sorted, scalarInsertionIndex = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a single-element scalar `insertion` with `index = [1,2]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         setInserted_D0_SK_ENABLED
            array = "AA"
            insertion = "XY"
            arrayNewS_ref = "XAXAX"
            arrayNewV_ref = "XYAXYAXY"
#elif       setInserted_D1_SK_ENABLED
            array = ["AA", "AA"]
            insertion = ["XX", "YY"]
            arrayNewS_ref = ["XX", "AA", "XX", "AA", "XX"]
            arrayNewV_ref = ["XX", "YY", "AA", "XX", "YY", "AA", "XX", "YY"]
#elif       setInserted_D1_IK_ENABLED
            array = [1_IKG, 1_IKG]
            insertion = [2_IKG, 3_IKG]
            arrayNewS_ref = [2_IKG, 1_IKG, 2_IKG, 1_IKG, 2_IKG]
            arrayNewV_ref = [2_IKG, 3_IKG, 1_IKG, 2_IKG, 3_IKG, 1_IKG, 2_IKG, 3_IKG]
#elif       setInserted_D1_CK_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            insertion = [(2._CKG,-2._CKG), (3._CKG,-3._CKG)]
            arrayNewS_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG)]
            arrayNewV_ref = [(2._CKG,-2._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG)]
#elif       setInserted_D1_RK_ENABLED
            array = [1._RKG, 1._RKG]
            insertion = [2._RKG, 3._RKG]
            arrayNewS_ref = [2._RKG, 1._RKG, 2._RKG, 1._RKG, 2._RKG]
            arrayNewV_ref = [2._RKG, 3._RKG, 1._RKG, 2._RKG, 3._RKG, 1._RKG, 2._RKG, 3._RKG]
#elif       setInserted_D1_LK_ENABLED
            array = [.false._LK, .false._LK]
            insertion = [.true._LK, .true._LK]
            arrayNewS_ref = [.true._LK, .false._LK, .true._LK, .false._LK, .true._LK]
            arrayNewV_ref = [.true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .true._LK, .true._LK]
#endif
            index = [1_IK, 2_IK, 3_IK]
            call runTestWith(positive, sorted)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a double-element vector `insertion` with `index = [1,2,3]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            call runTestWith(positive, sorted, scalarInsertionIndex = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a single-element scalar `insertion` with `index = [1,2,3]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (.not. getOption(.false._LK,positive)) then

                call reset()

#if             setInserted_D0_SK_ENABLED
                array = "AA"
                insertion = "XY"
                arrayNewS_ref = "XAXAXX"
                arrayNewV_ref = "XYAXYAXYXY"
#elif           setInserted_D1_SK_ENABLED
                array = ["AA", "AA"]
                insertion = ["XX", "YY"]
                arrayNewS_ref = ["XX", "AA", "XX", "AA", "XX", "XX"]
                arrayNewV_ref = ["XX", "YY", "AA", "XX", "YY", "AA", "XX", "YY", "XX", "YY"]
#elif           setInserted_D1_IK_ENABLED
                array = [1_IKG, 1_IKG]
                insertion = [2_IKG, 3_IKG]
                arrayNewS_ref = [2_IKG, 1_IKG, 2_IKG, 1_IKG, 2_IKG, 2_IKG]
                arrayNewV_ref = [2_IKG, 3_IKG, 1_IKG, 2_IKG, 3_IKG, 1_IKG, 2_IKG, 3_IKG, 2_IKG, 3_IKG]
#elif           setInserted_D1_CK_ENABLED
                array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
                insertion = [(2._CKG,-2._CKG), (3._CKG,-3._CKG)]
                arrayNewS_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (2._CKG,-2._CKG)]
                arrayNewV_ref = [(2._CKG,-2._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG)]
#elif           setInserted_D1_RK_ENABLED
                array = [1._RKG, 1._RKG]
                insertion = [2._RKG, 3._RKG]
                arrayNewS_ref = [2._RKG, 1._RKG, 2._RKG, 1._RKG, 2._RKG, 2._RKG]
                arrayNewV_ref = [2._RKG, 3._RKG, 1._RKG, 2._RKG, 3._RKG, 1._RKG, 2._RKG, 3._RKG, 2._RKG, 3._RKG]
#elif           setInserted_D1_LK_ENABLED
                array = [.false._LK, .false._LK]
                insertion = [.true._LK, .true._LK]
                arrayNewS_ref = [.true._LK, .false._LK, .true._LK, .false._LK, .true._LK, .true._LK]
                arrayNewV_ref = [.true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .true._LK, .true._LK]
#endif
                index = [1_IK, 2_IK, 3_IK, 0_IK]
                call runTestWith(positive, sorted)
                call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a double-element vector `insertion` with `index = [1_IK, 2_IK, 3_IK, 0_IK]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

                call runTestWith(positive, sorted, scalarInsertionIndex = 1_IK)
                call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a single-element scalar `insertion` with `index = [1_IK, 2_IK, 3_IK, 0_IK]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

                index = [-2_IK, -1_IK, 0_IK, 0_IK]
                call runTestWith(positive, sorted)
                call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a double-element vector `insertion` with `index = [-2_IK, -1_IK, 0_IK, 0_IK]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

                call runTestWith(positive, sorted, scalarInsertionIndex = 1_IK)
                call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a single-element scalar `insertion` with `index = [-2_IK, -1_IK, 0_IK, 0_IK]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (.not. getOption(.false._LK, sorted)) then

                call reset()

#if             setInserted_D0_SK_ENABLED
                array = "AA"
                insertion = "XY"
                arrayNewS_ref = "XAXAXX"
                arrayNewV_ref = "XYAXYAXYXY"
#elif           setInserted_D1_SK_ENABLED
                array = ["AA", "AA"]
                insertion = ["XX", "YY"]
                arrayNewS_ref = ["XX", "AA", "XX", "AA", "XX", "XX"]
                arrayNewV_ref = ["XX", "YY", "AA", "XX", "YY", "AA", "XX", "YY", "XX", "YY"]
#elif           setInserted_D1_IK_ENABLED
                array = [1_IKG, 1_IKG]
                insertion = [2_IKG, 3_IKG]
                arrayNewS_ref = [2_IKG, 1_IKG, 2_IKG, 1_IKG, 2_IKG, 2_IKG]
                arrayNewV_ref = [2_IKG, 3_IKG, 1_IKG, 2_IKG, 3_IKG, 1_IKG, 2_IKG, 3_IKG, 2_IKG, 3_IKG]
#elif           setInserted_D1_CK_ENABLED
                array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
                insertion = [(2._CKG,-2._CKG), (3._CKG,-3._CKG)]
                arrayNewS_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (2._CKG,-2._CKG)]
                arrayNewV_ref = [(2._CKG,-2._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG)]
#elif           setInserted_D1_RK_ENABLED
                array = [1._RKG, 1._RKG]
                insertion = [2._RKG, 3._RKG]
                arrayNewS_ref = [2._RKG, 1._RKG, 2._RKG, 1._RKG, 2._RKG, 2._RKG]
                arrayNewV_ref = [2._RKG, 3._RKG, 1._RKG, 2._RKG, 3._RKG, 1._RKG, 2._RKG, 3._RKG, 2._RKG, 3._RKG]
#elif           setInserted_D1_LK_ENABLED
                array = [.false._LK, .false._LK]
                insertion = [.true._LK, .true._LK]
                arrayNewS_ref = [.true._LK, .false._LK, .true._LK, .false._LK, .true._LK, .true._LK]
                arrayNewV_ref = [.true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .true._LK, .true._LK]
#endif
                index = [3_IK, 3_IK, 1_IK, 2_IK]
                call runTestWith(positive, sorted)
                call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a double-element vector `insertion` with `index = [3_IK, 0_IK, 2_IK, 1_IK]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

                call runTestWith(positive, sorted, scalarInsertionIndex = 1_IK)
                call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a single-element scalar `insertion` with `index = [3_IK, 0_IK, 2_IK, 1_IK]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

                if (.not. getOption(.false._LK, positive)) then
                    index = [-1_IK, 0_IK, -2_IK, 0_IK]
                    call runTestWith(positive, sorted)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a double-element vector `insertion` with `index = [-1_IK, 0_IK, -2_IK, 0_IK]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))

                    call runTestWith(positive, sorted, scalarInsertionIndex = 1_IK)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a single-element scalar `insertion` with `index = [-1_IK, 0_IK, -2_IK, 0_IK]` must yield a proper `arrayNew` with getInsertedEnabled = "//getStr(getInsertedEnabled)//".", int(__LINE__, IK))
                end if

            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(index)) deallocate(index)
            if (allocated(array)) deallocate(array)
            if (allocated(insertion)) deallocate(insertion)
            if (allocated(arrayNewV_ref)) deallocate(arrayNewV_ref)
        end subroutine reset

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestWith(positive, sorted, scalarInsertionIndex)
            use pm_arrayResize, only: setResized
            logical(LK) , intent(in), optional  :: positive, sorted
            integer(IK) , intent(in), optional  :: scalarInsertionIndex

            integer(IK) :: lenArrayNew
            lenArrayNew = GET_SIZE(array, kind = IK) + size(index, kind = IK)

            if (present(scalarInsertionIndex)) then
                call setResized(arrayNew, lenArrayNew)
                if (getInsertedEnabled) then
                    call setInserted(arrayNew, array, insertion(GET_INDEX(scalarInsertionIndex)), index = index, positive = positive, sorted = sorted)
                else
                    arrayNew = getInserted(array, insertion(GET_INDEX(scalarInsertionIndex)), index = index, positive = positive, sorted = sorted)
                end if
                assertion = assertion .and. ALL(arrayNew IS_EQUAL arrayNewS_ref)
                call reportFailure(positive, sorted)
            else
                lenArrayNew = lenArrayNew + size(index, kind = IK) * (GET_SIZE(insertion, kind = IK) - 1_IK)
                call setResized(arrayNew, lenArrayNew)
                if (getInsertedEnabled) then
                    call setInserted(arrayNew, array, insertion, index = index, positive = positive, sorted = sorted)
                else
                    arrayNew = getInserted(array, insertion, index = index, positive = positive, sorted = sorted)
                end if
                assertion = assertion .and. ALL(arrayNew IS_EQUAL arrayNewV_ref)
                call reportFailure(positive, sorted)
            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reportFailure(positive, sorted)

            use pm_io, only: display_type

            logical(LK) , intent(in), optional  :: positive, sorted

            type(display_type) :: disp
            disp = display_type()

            if (test%traceable .and. .not. assertion) then

                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")

                call disp%show("arrayNew")
                call disp%show( arrayNew )
                call disp%show("arrayNewV_ref")
                call disp%show( arrayNewV_ref )
                call disp%show("index")
                call disp%show( index )
                call disp%show("present(positive)")
                call disp%show( present(positive) )
                call disp%show("present(sorted)")
                call disp%show( present(sorted) )
    
                if (present(sorted)) then
                call disp%show("sorted")
                call disp%show( sorted )
                end if

                if (present(positive)) then
                call disp%show("positive")
                call disp%show( positive )
                end if

                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP

            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_EQUAL
#undef  ALL
