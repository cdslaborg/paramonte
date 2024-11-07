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
!>  [setRemoved](@ref pm_arrayRemove::setRemoved).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif
        character(*, SK), parameter :: PROCEDURE_NAME = "@setRemoved()"

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     setRemoved_ENABLED && D1_D0_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D1_D0_ENABLED
        character(2,SKG), dimension(:), allocatable :: array, arrayNew, arrayNew_ref
        character(2,SKG)                            :: pattern
#elif   IK_ENABLED && D1_D0_ENABLED
        integer(IKG)    , dimension(:), allocatable :: array, arrayNew, arrayNew_ref
        integer(IKG)                                :: pattern
#elif   LK_ENABLED && D1_D0_ENABLED
        logical(LKG)    , dimension(:), allocatable :: array, arrayNew, arrayNew_ref
        logical(LKG)                                :: pattern
#elif   CK_ENABLED && D1_D0_ENABLED
        complex(CKG)    , dimension(:), allocatable :: array, arrayNew, arrayNew_ref
        complex(CKG)                                :: pattern
#elif   RK_ENABLED && D1_D0_ENABLED
        real(RKG)       , dimension(:), allocatable :: array, arrayNew, arrayNew_ref
        real(RKG)                                   :: pattern
#else
#error  "Unrecognized interface."
#endif
        integer(IK) , allocatable   :: instance(:)
        logical(LK)                 :: getRemovedEnabled

        assertion = .true._LK

        getRemovedEnabled = .false._LK
        call runTestsWith()
        call runTestsWith(iseq = iseq)

        getRemovedEnabled = .true._LK
        call runTestsWith()
        call runTestsWith(iseq = iseq)

    contains

        function iseq(segment, pattern) result(equivalent)
#if         SK_ENABLED && D1_D0_ENABLED
            character(*,SKG), intent(in) :: segment, pattern
#elif       IK_ENABLED && D1_D0_ENABLED
            integer(IKG)    , intent(in) :: segment, pattern
#elif       CK_ENABLED && D1_D0_ENABLED
            complex(CKG)    , intent(in) :: segment, pattern
#elif       RK_ENABLED && D1_D0_ENABLED
            real(RKG)       , intent(in) :: segment, pattern
#elif       LK_ENABLED && D1_D0_ENABLED
            logical(LKG)    , intent(in) :: segment, pattern
#endif
            logical(LK) :: equivalent
            equivalent = pattern IS_EQUAL segment
        end function

        subroutine reset()
            if (allocated(array)) deallocate(array)
            if (allocated(instance)) deallocate(instance)
            if (allocated(arrayNew_ref)) deallocate(arrayNew_ref)
        end subroutine reset

        subroutine runTestsWith(iseq)
            logical(LK), external, optional  :: iseq

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            pattern = " "
#elif       IK_ENABLED && D1_D0_ENABLED
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            pattern = 1._CKG
#elif       RK_ENABLED && D1_D0_ENABLED
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            pattern = .false._LKG
#endif
            allocate(arrayNew_ref(0))
            allocate(instance(0))
            allocate(array(0))

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [1_IK, 2_IK]
            allocate(arrayNew_ref(0))

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .false._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [0_IK, 1_IK, 2_IK, 3_IK]
            allocate(arrayNew_ref(0))

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .false._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [1_IK, -1_IK]
            allocate(arrayNew_ref(0))

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .false._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [1_IK, -2_IK, -1_IK]
            allocate(arrayNew_ref(0))

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [1_IK, -1_IK, -2_IK, 0_IK]
            allocate(arrayNew_ref(0))

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, unsorted, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "BB", "AA"]
            pattern = "AA"
            arrayNew_ref = ["BB", "AA"]
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 2_IKG, 1_IKG]
            pattern = 1_IKG
            arrayNew_ref = [2_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
            arrayNew_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 2._RKG, 1._RKG]
            pattern = 1._RKG
            arrayNew_ref = [2._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .true._LK, .false._LK]
            pattern = .false._LK
            arrayNew_ref = [.true._LK, .false._LK]
#endif
            instance = [1_IK, -2_IK]

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "BB", "AA"]
            pattern = "AA"
            arrayNew_ref = ["BB", "AA"]
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 2_IKG, 1_IKG]
            pattern = 1_IKG
            arrayNew_ref = [2_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
            arrayNew_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 2._RKG, 1._RKG]
            pattern = 1._RKG
            arrayNew_ref = [2._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .true._LK, .false._LK]
            pattern = .false._LK
            arrayNew_ref = [.true._LK, .false._LK]
#endif
            instance = [-2_IK]

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the scalar `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["BB", "AA", "AA"]
            pattern = "AA"
            arrayNew_ref = ["BB", "AA"]
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [2_IKG, 1_IKG, 1_IKG]
            pattern = 1_IKG
            arrayNew_ref = [2_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(2._CKG,-2._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
            arrayNew_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [2._RKG, 1._RKG, 1._RKG]
            pattern = 1._RKG
            arrayNew_ref = [2._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.true._LK, .false._LK, .false._LK]
            pattern = .false._LK
            arrayNew_ref = [.true._LK, .false._LK]
#endif
            instance = [-1_IK, 3_IK]

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the scalar `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the scalar `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the scalar `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the scalar `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the scalar `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the scalar `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the scalar `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the scalar `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA"]
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG]
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG]
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK]
#endif
            allocate(arrayNew_ref(0))
            pattern = array(1)
            instance = [-1_IK, 3_IK]

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA"]
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG]
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG]
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.true._LKG] ! \warning \todo The GNU 10-11 fail on this test if array = [.false._LKG] for LK5_ENABLED. This is super odd and needs future investigation. ifort does not have LK5 but gracefully passes all other tests.
#endif
            pattern = array(1)
            instance = [0_IK, 3_IK]
            arrayNew_ref = array

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            allocate(character(2) :: array(0))
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            allocate(array(0))
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            allocate(array(0))
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            allocate(array(0))
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            allocate(array(0))
            pattern = .false._LK
#endif
            arrayNew_ref = array
            instance = [0_IK, 3_IK]

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "BB", "CC", "AA"]
            pattern = "AA"
            arrayNew_ref = ["BB", "CC"]
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 2_IKG, 2_IKG, 3_IKG, 3_IKG, 1_IKG]
            pattern = 1_IKG
            arrayNew_ref = [2_IKG, 2_IKG, 3_IKG, 3_IKG]
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
            arrayNew_ref = [(2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG)]
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 2._RKG, 2._RKG, 3._RKG, 3._RKG, 1._RKG]
            pattern = 1._RKG
            arrayNew_ref = [2._RKG, 2._RKG, 3._RKG, 3._RKG]
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .true._LK, .true._LK, .true._LK, .true._LK, .false._LK]
            pattern = .false._LK
            arrayNew_ref = [.true._LK, .true._LK, .true._LK, .true._LK]
#endif
            instance = [1_IK, -1_IK]

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            pattern = "XX"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = 0_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (0._CKG,0._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = 0._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            pattern = .true._LK
#endif
            instance = [0_IK, 3_IK]
            arrayNew_ref = array

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "XX", "AA", "XX", "AA"]
            pattern = "XX"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 0_IKG, 1_IKG, 0_IKG, 1_IKG]
            pattern = 0_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG)]
            pattern = (0._CKG,0._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 0._RKG, 1._RKG, 0._RKG, 1._RKG]
            pattern = 0._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK]
            pattern = .true._LK
#endif
            instance = [1_IK, -1_IK]
            allocate(arrayNew_ref(3))
            arrayNew_ref(1:1) = array(1:1)
            arrayNew_ref(2:2) = array(3:3)
            arrayNew_ref(3:3) = array(5:5)

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [0_IK, 3_IK]
            arrayNew_ref = array

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        subroutine report   ( iseq & ! LCOV_EXCL_LINE
                            , instance & ! LCOV_EXCL_LINE
                            , sorted & ! LCOV_EXCL_LINE
                            , unique & ! LCOV_EXCL_LINE
                            )
            logical(LK) , external  , optional              :: iseq
            integer(IK) , intent(in), optional, contiguous  :: instance(:)
            logical(LK) , intent(in), optional              :: sorted
            logical(LK) , intent(in), optional              :: unique

            if (getRemovedEnabled) then
                arrayNew = array
                if (present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq, instance = instance, sorted = sorted, unique = unique)
                elseif (present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq, instance = instance, sorted = sorted)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq, instance = instance, unique = unique)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq, instance = instance)
                elseif (present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                    call setRemoved(arrayNew, pattern, instance = instance, sorted = sorted, unique = unique)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, instance = instance, sorted = sorted)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                    call setRemoved(arrayNew, pattern, instance = instance, unique = unique)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, instance = instance)
                elseif (.not. present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern)
                else
                    error stop PROCEDURE_NAME//SK_": Unrecognized interface in testing." ! LCOV_EXCL_LINE
                end if
            else
                if (present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                    arrayNew = getRemoved(array, pattern, iseq = iseq, instance = instance, sorted = sorted, unique = unique)
                elseif (present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(array, pattern, iseq = iseq, instance = instance, sorted = sorted)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                    arrayNew = getRemoved(array, pattern, iseq = iseq, instance = instance, unique = unique)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(array, pattern, iseq = iseq, instance = instance)
                elseif (present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(array, pattern, iseq = iseq)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                    arrayNew = getRemoved(array, pattern, instance = instance, sorted = sorted, unique = unique)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(array, pattern, instance = instance, sorted = sorted)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                    arrayNew = getRemoved(array, pattern, instance = instance, unique = unique)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(array, pattern, instance = instance)
                elseif (.not. present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(array, pattern)
                else
                    error stop PROCEDURE_NAME//SK_": Unrecognized interface in testing." ! LCOV_EXCL_LINE
                end if
            end if

            ! Report test results if needed.

            assertion = assertion .and. all(arrayNew IS_EQUAL arrayNew_ref)

            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("arrayNew")
                call test%disp%show( arrayNew )
                call test%disp%show("arrayNew_ref")
                call test%disp%show( arrayNew_ref )
                call test%disp%show("present(iseq)")
                call test%disp%show( present(iseq) )
                call test%disp%show("present(instance)")
                call test%disp%show( present(instance) )
                if (present(instance)) then
                    call test%disp%show("instance")
                    call test%disp%show( instance )
                end if
                call test%disp%show("present(unique)")
                call test%disp%show( present(unique) )
                if (present(unique)) then
                    call test%disp%show("unique")
                    call test%disp%show( unique )
                end if
                call test%disp%show("present(sorted)")
                call test%disp%show( present(sorted) )
                if (present(sorted)) then
                    call test%disp%show("sorted")
                    call test%disp%show( sorted )
                end if
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRemoved_ENABLED && (D0_D0_ENABLED || D1_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_D0_ENABLED
        character(:,SKG), allocatable               :: Array, pattern, arrayNew, arrayNew_ref
#elif   SK_ENABLED && D1_D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: Array, pattern, arrayNew, arrayNew_ref
#elif   IK_ENABLED && D1_D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: Array, pattern, arrayNew, arrayNew_ref
#elif   LK_ENABLED && D1_D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: Array, pattern, arrayNew, arrayNew_ref
#elif   CK_ENABLED && D1_D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: Array, pattern, arrayNew, arrayNew_ref
#elif   RK_ENABLED && D1_D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: Array, pattern, arrayNew, arrayNew_ref
#else
#error  "Unrecognized interface."
#endif
        integer(IK) , allocatable   :: instance(:)
        logical(LK)                 :: getRemovedEnabled

        assertion = .true._LK

        getRemovedEnabled = .false._LK
        call runTestsWith()
        call runTestsWith(iseq = iseq)

        getRemovedEnabled = .true._LK
        call runTestsWith()
        call runTestsWith(iseq = iseq)

    contains

        function iseq(segment, pattern & ! LCOV_EXCL_LINE
#if     D1_D1_ENABLED
        , lenPattern & ! LCOV_EXCL_LINE
#endif
        ) result(equivalent)
#if     D1_D1_ENABLED
            integer(IK), intent(in) :: lenPattern
#endif
#if         SK_ENABLED && D0_D0_ENABLED
            character(*,SKG), intent(in) :: segment, pattern
#elif       SK_ENABLED && D1_D1_ENABLED
            character(*,SKG), intent(in) :: segment(lenPattern), pattern(lenPattern)
#elif       IK_ENABLED && D1_D1_ENABLED
            integer(IKG)    , intent(in) :: segment(lenPattern), pattern(lenPattern)
#elif       CK_ENABLED && D1_D1_ENABLED
            complex(CKG)    , intent(in) :: segment(lenPattern), pattern(lenPattern)
#elif       RK_ENABLED && D1_D1_ENABLED
            real(RKG)       , intent(in) :: segment(lenPattern), pattern(lenPattern)
#elif       LK_ENABLED && D1_D1_ENABLED
            logical(LKG)    , intent(in) :: segment(lenPattern), pattern(lenPattern)
#endif
            logical(LK)                  :: equivalent
            equivalent = all([pattern IS_EQUAL segment])
        end function

        subroutine reset()
            if (allocated(Array)) deallocate(Array)
            if (allocated(pattern)) deallocate(pattern)
            if (allocated(instance)) deallocate(instance)
            if (allocated(arrayNew_ref)) deallocate(arrayNew_ref)
        end subroutine reset

        subroutine runTestsWith(iseq)
            logical(LK), external, optional  :: iseq

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            pattern = " "
            allocate(character(0,SKG) :: arrayNew_ref, Array)
#elif       SK_ENABLED && D1_D1_ENABLED
            pattern = [" "]
            allocate(character(2,SKG) :: arrayNew_ref(0), Array(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            pattern = [1_IKG]
            allocate(arrayNew_ref(0), Array(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            pattern = [1._CKG]
            allocate(arrayNew_ref(0), Array(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            pattern = [1._RKG]
            allocate(arrayNew_ref(0), Array(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            pattern = [.false._LK]
            allocate(arrayNew_ref(0), Array(0))
#endif
            allocate(instance(0))

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AAAA"
            pattern = "AA"
            allocate(character(0,SKG) :: arrayNew_ref)
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            pattern = ["AA"]
            allocate(character(2,SKG) :: arrayNew_ref(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
            allocate(arrayNew_ref(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
            allocate(arrayNew_ref(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
            allocate(arrayNew_ref(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            pattern = [.false._LK]
            allocate(arrayNew_ref(0))
#endif
            instance = [1_IK, 2_IK]

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .false._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AAAA"
            pattern = "AA"
            allocate(character(0,SKG) :: arrayNew_ref)
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            pattern = ["AA"]
            allocate(character(2,SKG) :: arrayNew_ref(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
            allocate(arrayNew_ref(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
            allocate(arrayNew_ref(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
            allocate(arrayNew_ref(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            pattern = [.false._LK]
            allocate(arrayNew_ref(0))
#endif
            instance = [0_IK, 1_IK, 2_IK, 3_IK]

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .false._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and redundant with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AAAA"
            pattern = "AA"
            allocate(character(0,SKG) :: arrayNew_ref)
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            pattern = ["AA"]
            allocate(character(2,SKG) :: arrayNew_ref(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
            allocate(arrayNew_ref(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
            allocate(arrayNew_ref(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
            allocate(arrayNew_ref(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            pattern = [.false._LK]
            allocate(arrayNew_ref(0))
#endif
            instance = [1_IK, -1_IK]

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .false._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty and mixed with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AAAA"
            pattern = "AA"
            allocate(character(0,SKG) :: arrayNew_ref)
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            pattern = ["AA"]
            allocate(character(2,SKG) :: arrayNew_ref(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
            allocate(arrayNew_ref(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
            allocate(arrayNew_ref(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
            allocate(arrayNew_ref(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            pattern = [.false._LK]
            allocate(arrayNew_ref(0))
#endif
            instance = [1_IK, -2_IK, -1_IK]

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AAAA"
            pattern = "AA"
            allocate(character(0,SKG) :: arrayNew_ref)
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            pattern = ["AA"]
            allocate(character(2,SKG) :: arrayNew_ref(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
            allocate(arrayNew_ref(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
            allocate(arrayNew_ref(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
            allocate(arrayNew_ref(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            pattern = [.false._LK]
            allocate(arrayNew_ref(0))
#endif
            instance = [1_IK, -1_IK, -2_IK, 0_IK]

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, unsorted, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield an empty `arrayNew` when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AABBAA"
            pattern = "AA"
            arrayNew_ref = "BBAA"
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "BB", "AA"]
            pattern = ["AA"]
            arrayNew_ref = ["BB", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 2_IKG, 1_IKG]
            pattern = [1_IKG]
            arrayNew_ref = [2_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
            arrayNew_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 2._RKG, 1._RKG]
            pattern = [1._RKG]
            arrayNew_ref = [2._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .true._LK, .false._LK]
            pattern = [.false._LK]
            arrayNew_ref = [.true._LK, .false._LK]
#endif
            instance = [1_IK, -2_IK]

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AABBAA"
            pattern = "AA"
            arrayNew_ref = "BBAA"
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "BB", "AA"]
            pattern = ["AA"]
            arrayNew_ref = ["BB", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 2_IKG, 1_IKG]
            pattern = [1_IKG]
            arrayNew_ref = [2_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
            arrayNew_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 2._RKG, 1._RKG]
            pattern = [1._RKG]
            arrayNew_ref = [2._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .true._LK, .false._LK]
            pattern = [.false._LK]
            arrayNew_ref = [.true._LK, .false._LK]
#endif
            instance = [-2_IK]

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the beginning and the end matching the single-element vector `pattern` must yield an `arrayNew` with the last two elements of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "BBAAAA"
            pattern = "AA"
            arrayNew_ref = "BBAA"
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["BB", "AA", "AA"]
            pattern = ["AA"]
            arrayNew_ref = ["BB", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [2_IKG, 1_IKG, 1_IKG]
            pattern = [1_IKG]
            arrayNew_ref = [2_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(2._CKG,-2._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
            arrayNew_ref = [(2._CKG,-2._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [2._RKG, 1._RKG, 1._RKG]
            pattern = [1._RKG]
            arrayNew_ref = [2._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.true._LK, .false._LK, .false._LK]
            pattern = [.false._LK]
            arrayNew_ref = [.true._LK, .false._LK]
#endif
            instance = [-1_IK, 3_IK]

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield an `arrayNew` of size two of the first two elements of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AA"
            allocate(character(0,SKG) :: arrayNew_ref)
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA"]
            allocate(character(2,SKG) :: arrayNew_ref(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG]
            allocate(arrayNew_ref(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG)]
            allocate(arrayNew_ref(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG]
            allocate(arrayNew_ref(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK]
            allocate(arrayNew_ref(0))
#endif
            pattern = Array
            instance = [-1_IK, 3_IK]

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK]
#endif
            pattern = Array
            instance = [0_IK, 3_IK]
            arrayNew_ref = Array

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            allocate(character(0,SKG) :: Array)
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            allocate(character(2,SKG) :: Array(0))
            pattern = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            allocate(Array(0))
            pattern = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            allocate(Array(0))
            pattern = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            allocate(Array(0))
            pattern = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            allocate(Array(0))
            pattern = [.false._LK]
#endif
            arrayNew_ref = Array
            instance = [0_IK, 3_IK]

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as the `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AABBCCAA"
            pattern = "AA"
            arrayNew_ref = "BBCC"
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "BB", "CC", "AA"]
            pattern = ["AA"]
            arrayNew_ref = ["BB", "CC"]
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 2_IKG, 2_IKG, 3_IKG, 3_IKG, 1_IKG]
            pattern = [1_IKG]
            arrayNew_ref = [2_IKG, 2_IKG, 3_IKG, 3_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
            arrayNew_ref = [(2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 2._RKG, 2._RKG, 3._RKG, 3._RKG, 1._RKG]
            pattern = [1._RKG]
            arrayNew_ref = [2._RKG, 2._RKG, 3._RKG, 3._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .true._LK, .true._LK, .true._LK, .true._LK, .false._LK]
            pattern = [.false._LK]
            arrayNew_ref = [.true._LK, .true._LK, .true._LK, .true._LK]
#endif
            instance = [1_IK, -1_IK]

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `arrayNew` of the remaining elements of `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AAAA"
            pattern = "XX"
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            pattern = ["XX"]
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            pattern = [0_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(0._CKG,0._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            pattern = [0._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            pattern = [.true._LK]
#endif
            instance = [0_IK, 3_IK]
            arrayNew_ref = Array

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AXAXA"
            pattern = "X"
            allocate(character(3) :: arrayNew_ref)
#elif       SK_ENABLED && D1_D1_ENABLED
            allocate(arrayNew_ref(3))
            Array = ["AA", "XX", "AA", "XX", "AA"]
            pattern = ["XX"]
#elif       IK_ENABLED && D1_D1_ENABLED
            allocate(arrayNew_ref(3))
            Array = [1_IKG, 0_IKG, 1_IKG, 0_IKG, 1_IKG]
            pattern = [0_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            allocate(arrayNew_ref(3))
            Array = [(1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG)]
            pattern = [(0._CKG,0._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            allocate(arrayNew_ref(3))
            Array = [1._RKG, 0._RKG, 1._RKG, 0._RKG, 1._RKG]
            pattern = [0._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            allocate(arrayNew_ref(3))
            Array = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK]
            pattern = [.true._LK]
#endif
            instance = [1_IK, -1_IK]
            arrayNew_ref(1:1) = Array(1:1)
            arrayNew_ref(2:2) = Array(3:3)
            arrayNew_ref(3:3) = Array(5:5)

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield an `arrayNew` of size three whose elements are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AAAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            pattern = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            pattern = [.false._LK]
#endif
            instance = [0_IK, 3_IK]
            arrayNew_ref = Array

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an `arrayNew` of size one whose value is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.` with getRemovedEnabled = "//getStr(getRemovedEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AA"
            pattern = ""
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            allocate(character(2,SKG) :: pattern(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            allocate(pattern(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [1._CKG, 1._CKG]
            allocate(pattern(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            allocate(pattern(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            allocate(pattern(0))
#endif
            arrayNew_ref = Array
            allocate(instance(0))

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `Array`.", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `array` when `instance` is empty.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `array` when `instance` is empty with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `array` when `instance` is empty with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `array` when `instance` is empty with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `array` when `instance` is empty with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `array` when `instance` is empty with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `array` when `instance` is empty with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `array` when `instance` is empty with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` that is the same as `array` when `instance` is empty with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AA"
            pattern = ""
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            allocate(character(2,SKG) :: pattern(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            allocate(pattern(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [1._CKG, 1._CKG]
            allocate(pattern(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            allocate(pattern(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            allocate(pattern(0))
#endif
            arrayNew_ref = Array
            instance = [1_IK]

            call report(iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `Array`.", int(__LINE__, IK))

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `array` when `instance` is non-empty.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `array` when `instance` is non-empty with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `array` when `instance` is non-empty with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `array` when `instance` is non-empty with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `array` when `instance` is non-empty with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `array` when `instance` is non-empty with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `array` when `instance` is non-empty with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `array` when `instance` is non-empty with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield an `arrayNew` whose value is the same as `array` when `instance` is non-empty with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AA"
            allocate(character(0,SKG) :: arrayNew_ref)
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            allocate(arrayNew_ref(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            allocate(arrayNew_ref(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            allocate(arrayNew_ref(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            allocate(arrayNew_ref(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            allocate(arrayNew_ref(0))
#endif
            pattern = Array
            instance = [-1_IK, 3_IK]

            call report(iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `arrayNew` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
#endif
            pattern = Array
            instance = [0_IK, 3_IK]
            arrayNew_ref = Array

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            Array = "AA"
            pattern = Array//Array
#elif       SK_ENABLED && D1_D1_ENABLED
            Array = ["AA", "AA"]
            pattern = [Array, Array]
#elif       IK_ENABLED && D1_D1_ENABLED
            Array = [1_IKG, 1_IKG]
            pattern = [Array, Array]
#elif       CK_ENABLED && D1_D1_ENABLED
            Array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [Array, Array]
#elif       RK_ENABLED && D1_D1_ENABLED
            Array = [1._RKG, 1._RKG]
            pattern = [Array, Array]
#elif       LK_ENABLED && D1_D1_ENABLED
            Array = [.false._LK, .false._LK]
            pattern = [Array, Array]
#endif
            instance = [0_IK, 3_IK]
            arrayNew_ref = Array

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            instance = [1_IK]

            call report(iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an `arrayNew` that is the same as `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        subroutine report   ( iseq & ! LCOV_EXCL_LINE
                            , instance & ! LCOV_EXCL_LINE
                            , sorted & ! LCOV_EXCL_LINE
                            , unique & ! LCOV_EXCL_LINE
                            )
            logical(LK) , external  , optional              :: iseq
            integer(IK) , intent(in), optional, contiguous  :: instance(:)
            logical(LK) , intent(in), optional              :: sorted
            logical(LK) , intent(in), optional              :: unique

            if (getRemovedEnabled) then
                arrayNew = Array
                if (present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq, instance = instance, sorted = sorted, unique = unique)
                elseif (present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq, instance = instance, sorted = sorted)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq, instance = instance, unique = unique)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq, instance = instance)
                elseif (present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, iseq = iseq)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                    call setRemoved(arrayNew, pattern, instance = instance, sorted = sorted, unique = unique)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, instance = instance, sorted = sorted)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                    call setRemoved(arrayNew, pattern, instance = instance, unique = unique)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern, instance = instance)
                elseif (.not. present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    call setRemoved(arrayNew, pattern)
                else
                    error stop PROCEDURE_NAME//SK_": Unrecognized interface in testing." ! LCOV_EXCL_LINE
                end if
            else
                if (present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                    arrayNew = getRemoved(Array, pattern, iseq = iseq, instance = instance, sorted = sorted, unique = unique)
                elseif (present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(Array, pattern, iseq = iseq, instance = instance, sorted = sorted)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                    arrayNew = getRemoved(Array, pattern, iseq = iseq, instance = instance, unique = unique)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(Array, pattern, iseq = iseq, instance = instance)
                elseif (present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(Array, pattern, iseq = iseq)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                    arrayNew = getRemoved(Array, pattern, instance = instance, sorted = sorted, unique = unique)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(Array, pattern, instance = instance, sorted = sorted)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                    arrayNew = getRemoved(Array, pattern, instance = instance, unique = unique)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(Array, pattern, instance = instance)
                elseif (.not. present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                    arrayNew = getRemoved(Array, pattern)
                else
                    error stop PROCEDURE_NAME//SK_": Unrecognized interface in testing." ! LCOV_EXCL_LINE
                end if
            end if

            ! Report test results if needed.

            assertion = assertion .and. all([arrayNew IS_EQUAL arrayNew_ref])

            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("arrayNew")
                call test%disp%show( arrayNew )
                call test%disp%show("arrayNew_ref")
                call test%disp%show( arrayNew_ref )
                call test%disp%show("present(iseq)")
                call test%disp%show( present(iseq) )
                call test%disp%show("present(instance)")
                call test%disp%show( present(instance) )
                if (present(instance)) then
                    call test%disp%show("instance")
                    call test%disp%show( instance )
                end if
                call test%disp%show("present(unique)")
                call test%disp%show( present(unique) )
                if (present(unique)) then
                    call test%disp%show("unique")
                    call test%disp%show( unique )
                end if
                call test%disp%show("present(sorted)")
                call test%disp%show( present(sorted) )
                if (present(sorted)) then
                    call test%disp%show("sorted")
                    call test%disp%show( sorted )
                end if
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if

        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef IS_EQUAL