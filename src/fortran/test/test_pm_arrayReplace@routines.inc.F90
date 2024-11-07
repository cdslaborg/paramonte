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
!>  This module contains implementations of the tests of the procedures under the generic interfaces
!>  [getReplaced](@ref Array_pmod::getReplaced),
!>  [setReplaced](@ref Array_pmod::setReplaced).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED || \
        getReplaced_D1_D1_D0_LK_ENABLED || setReplaced_D1_D1_D0_LK_ENABLED || \
        getReplaced_D1_D0_D1_LK_ENABLED || setReplaced_D1_D0_D1_LK_ENABLED || \
        getReplaced_D1_D0_D0_LK_ENABLED || setReplaced_D1_D0_D0_LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
#define GET_SIZE len
#else
#define GET_SIZE size
#endif

#if     getReplaced_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@getReplaced()"
#elif   setReplaced_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@setReplaced()"
#endif

        integer(IK), allocatable :: instance(:)

#if     getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
#define ALL
        character(:,SKG), allocatable   :: ArrayNew_ref, ArrayTemplate, Array, pattern, Replacement
#elif   getReplaced_D1_D0_D0_SK_ENABLED || setReplaced_D1_D0_D0_SK_ENABLED
        character(2,SKG), allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:)
        character(2,SKG)                :: pattern, Replacement, lower, upper
        lower = SKG_"aa"; upper = SKG_"zz"
#elif   getReplaced_D1_D0_D0_IK_ENABLED || setReplaced_D1_D0_D0_IK_ENABLED
        integer(IKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:)
        integer(IKG)                    :: pattern, Replacement, lower, upper
        lower = -127_IKG; upper = 127_IKG
#elif   getReplaced_D1_D0_D0_LK_ENABLED || setReplaced_D1_D0_D0_LK_ENABLED
        logical(LKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:)
        logical(LKG)                    :: pattern, Replacement, lower, upper
        lower = .false._LKG; upper = .true._LKG
#elif   getReplaced_D1_D0_D0_CK_ENABLED || setReplaced_D1_D0_D0_CK_ENABLED
        complex(CKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:)
        complex(CKG)                    :: pattern, Replacement, lower, upper
        lower = (0._CKG, 0._CKG); upper = (1._CKG, 1._CKG)
#elif   getReplaced_D1_D0_D0_RK_ENABLED || setReplaced_D1_D0_D0_RK_ENABLED
        real(RKG)       , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:)
        real(RKG)                       :: pattern, Replacement, lower, upper
        lower = 0._RKG; upper = 1._RKG
#elif   getReplaced_D1_D0_D1_SK_ENABLED || setReplaced_D1_D0_D1_SK_ENABLED
        character(2,SKG), allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), Replacement(:)
        character(2,SKG)                :: pattern, lower, upper
        lower = SKG_"aa"; upper = SKG_"zz"
#elif   getReplaced_D1_D0_D1_IK_ENABLED || setReplaced_D1_D0_D1_IK_ENABLED
        integer(IKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), Replacement(:)
        integer(IKG)                    :: pattern, lower, upper
        lower = 0_IKG; upper = 9_IKG
#elif   getReplaced_D1_D0_D1_LK_ENABLED || setReplaced_D1_D0_D1_LK_ENABLED
        logical(LKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), Replacement(:)
        logical(LKG)                    :: pattern, lower, upper
        lower = .false._LKG; upper = .true._LKG
#elif   getReplaced_D1_D0_D1_CK_ENABLED || setReplaced_D1_D0_D1_CK_ENABLED
        complex(CKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), Replacement(:)
        complex(CKG)                    :: pattern, lower, upper
        lower = (0._CKG, 0._CKG); upper = (1._CKG, 1._CKG)
#elif   getReplaced_D1_D0_D1_RK_ENABLED || setReplaced_D1_D0_D1_RK_ENABLED
        real(RKG)       , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), Replacement(:)
        real(RKG)                       :: pattern, lower, upper
        lower = 0._RKG; upper = 1._RKG
#elif   getReplaced_D1_D1_D0_SK_ENABLED || setReplaced_D1_D1_D0_SK_ENABLED
        character(2,SKG), allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:)
        character(2,SKG)                :: Replacement, lower, upper
        lower = SKG_"aa"; upper = SKG_"zz"
#elif   getReplaced_D1_D1_D0_IK_ENABLED || setReplaced_D1_D1_D0_IK_ENABLED
        integer(IKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:)
        integer(IKG)                    :: Replacement, lower, upper
        lower = 0_IKG; upper = 9_IKG
#elif   getReplaced_D1_D1_D0_LK_ENABLED || setReplaced_D1_D1_D0_LK_ENABLED
        logical(LKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:)
        logical(LKG)                    :: Replacement, lower, upper
        lower = .false._LKG; upper = .true._LKG
#elif   getReplaced_D1_D1_D0_CK_ENABLED || setReplaced_D1_D1_D0_CK_ENABLED
        complex(CKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:)
        complex(CKG)                    :: Replacement, lower, upper
        lower = (0._CKG, 0._CKG); upper = (1._CKG, 1._CKG)
#elif   getReplaced_D1_D1_D0_RK_ENABLED || setReplaced_D1_D1_D0_RK_ENABLED
        real(RKG)       , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:)
        real(RKG)                       :: Replacement, lower, upper
        lower = 0._RKG; upper = 1._RKG
#elif   getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
        character(2,SKG), allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:), Replacement(:)
#elif   getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
        integer(IKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:), Replacement(:)
#elif   getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
        logical(LKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:), Replacement(:)
#elif   getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
        complex(CKG)    , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:), Replacement(:)
#elif   getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
        real(RKG)       , allocatable   :: ArrayNew_ref(:), ArrayTemplate(:), Array(:), pattern(:), Replacement(:)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        call runTestsWith(getReplacedEnabled = .true._LK)
        call runTestsWith(getReplacedEnabled = .false._LK)
        call runTestsWith(getReplacedEnabled = .true._LK, iseq = iseq)
        call runTestsWith(getReplacedEnabled = .false._LK, iseq = iseq)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(getReplacedEnabled, iseq)
            logical(LK), intent(in)         :: getReplacedEnabled
            logical(LK), external, optional :: iseq

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if         getReplaced_D1_D0_D0_ENABLED || getReplaced_D1_D0_D1_ENABLED || getReplaced_D1_D1_D0_ENABLED || \
            setReplaced_D1_D0_D0_ENABLED || setReplaced_D1_D0_D1_ENABLED || setReplaced_D1_D1_D0_ENABLED
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            integer(IK) :: itest, jtest
            integer(IK), allocatable :: InstanceCopy(:)!, InstancePositive(:), Loc(:)

            do itest = 1_IK, 10_IK

                call reset()
                allocate(ArrayTemplate(getUnifRand(1, 10)))
                call setUnifRand(ArrayTemplate, lower, upper)
#if             getReplaced_D1_D1_D0_ENABLED || setReplaced_D1_D1_D0_ENABLED
                allocate(pattern(getUnifRand(0, size(ArrayTemplate) + 2)))
#endif
#if             getReplaced_D1_D0_D1_ENABLED || setReplaced_D1_D0_D1_ENABLED
                allocate(Replacement(getUnifRand(0, size(ArrayTemplate) + 2)))
#endif
                call setUnifRand(Replacement, lower, upper)

                do jtest = 1_IK, 2_IK

                    if (jtest == 1_IK) then
                        call setUnifRand(pattern, lower, upper)
                    else
#if                     getReplaced_D1_D1_D0_ENABLED || setReplaced_D1_D1_D0_ENABLED
                        pattern = ArrayTemplate(getUnifRand(1, size(ArrayTemplate) + 1) : getUnifRand(0, size(ArrayTemplate)))
#else
                        pattern = ArrayTemplate(getUnifRand(1, size(ArrayTemplate)))
#endif
                    end if

                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array.", int(__LINE__, IK))

                    if (allocated(instance)) deallocate(instance)
                    allocate(instance(getUnifRand(0, 2 * size(ArrayTemplate))))
                    call setUnifRand(instance, -2_IK * size(ArrayTemplate, kind = IK), 2_IK * size(ArrayTemplate, kind = IK))

                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq, instance = instance)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array with `instance`.", int(__LINE__, IK))

                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array with `instance` with `sorted = .false.`.", int(__LINE__, IK))

                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array with `instance` with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array with `instance` with `unique = .false.`.", int(__LINE__, IK))

                    !   Sorting and removing duplicates in mixed `instance` is very difficult without a priori knowledge of the locations of `pattern` in `Array`.
                    !   So, we limit the tests to only cases where `instance` is strictly positive or negative.

                    instance = abs(instance)
                    if (getUnifRand()) instance = -instance
                    InstanceCopy = instance
                    instance = getUnique(instance)

                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array with `instance` with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array with `instance` with `unique = .true.`.", int(__LINE__, IK))

                    call setSorted(InstanceCopy)
                    instance = InstanceCopy

                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array with `instance` with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array with `instance` with `sorted = .true.`.", int(__LINE__, IK))

                    instance = getUnique(instance)
                    Array = ArrayTemplate
                    call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": All instances of pattern must be properly replaced with the input replacement in the input array with `instance` with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

                end do

            end do


            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif       getReplaced_D0_D0_D0_ENABLED || getReplaced_D1_D1_D1_ENABLED || \
            setReplaced_D0_D0_D0_ENABLED || setReplaced_D1_D1_D1_ENABLED
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = ""
            pattern = ""
            Replacement = ""
            ArrayNew_ref = ""
#elif       getReplaced_D1_D1_D1_SK_ENABLED
            allocate(character(2,SKG) :: ArrayTemplate(0), ArrayNew_ref(0), Replacement(0), pattern(0))
#elif       setReplaced_D1_D1_D1_SK_ENABLED
            allocate(ArrayTemplate(0), ArrayNew_ref(0), Replacement(0), pattern(0))
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            allocate(ArrayTemplate(0), ArrayNew_ref(0), Replacement(0), pattern(0))
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            allocate(ArrayTemplate(0), ArrayNew_ref(0), Replacement(0), pattern(0))
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            allocate(ArrayTemplate(0), ArrayNew_ref(0), Replacement(0), pattern(0))
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            allocate(ArrayTemplate(0), ArrayNew_ref(0), Replacement(0), pattern(0))
#endif
            allocate(instance(0))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array with `instance`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array with `instance` with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array with `instance` with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array with `instance` with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array with `instance` with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array with `instance` with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array with `instance` with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array with `instance` with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, pattern, replacement, has empty resulting array with `instance` with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = ""
            pattern = ""
            Replacement = "XX"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            allocate(character(2,SKG) :: ArrayTemplate(0), pattern(0))
            Replacement = ["XX"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            allocate(ArrayTemplate(0), pattern(0))
            Replacement = [11_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            allocate(ArrayTemplate(0), pattern(0))
            Replacement = [11._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            allocate(ArrayTemplate(0), pattern(0))
            Replacement = [11._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            allocate(ArrayTemplate(0), pattern(0))
            Replacement = [.false._LKG]
#endif
            ArrayNew_ref = Replacement
            allocate(instance(0))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in replacement as new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            ArrayNew_ref = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance` is empty.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance` is empty with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance` is empty with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance` is empty with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance` is empty with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance` is empty with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance` is empty with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance` is empty with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance` is empty with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = ""
            pattern = ""
            Replacement = "XX"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            allocate(character(2,SKG) :: ArrayTemplate(0), pattern(0))
            Replacement = ["XX"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            allocate(ArrayTemplate(0), pattern(0))
            Replacement = [11_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            allocate(ArrayTemplate(0), pattern(0))
            Replacement = [11._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            allocate(ArrayTemplate(0), pattern(0))
            Replacement = [11._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            allocate(ArrayTemplate(0), pattern(0))
            Replacement = [.false._LKG]
#endif
            ArrayNew_ref = Replacement
            instance = [1000_IK]

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in replacement as new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            ArrayNew_ref = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance`is a nonsensical number.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance`is a nonsensical number with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance`is a nonsensical number with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance`is a nonsensical number with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance`is a nonsensical number with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance`is a nonsensical number with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance`is a nonsensical number with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance`is a nonsensical number with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with empty pattern and non-empty replacement, must result in an empty array when `instance`is a nonsensical number with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = "A"
            pattern = ""
            Replacement = "XX"
            ArrayNew_ref = "XXAXX"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            ArrayTemplate = ["AA"]
            allocate(character(2,SKG) :: pattern(0))
            Replacement = ["XX"]
            ArrayNew_ref = ["XX", "AA", "XX"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            ArrayTemplate = [33_IKG]
            allocate(pattern(0))
            Replacement = [11_IKG]
            ArrayNew_ref = [11_IKG, 33_IKG, 11_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            ArrayTemplate = [33._CKG]
            allocate(pattern(0))
            Replacement = [11._CKG]
            ArrayNew_ref = [11._CKG, 33._CKG, 11._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            ArrayTemplate = [33._RKG]
            allocate(pattern(0))
            Replacement = [11._RKG]
            ArrayNew_ref = [11._RKG, 33._RKG, 11._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            ArrayTemplate = [.true._LKG]
            allocate(pattern(0))
            Replacement = [.false._LKG]
            ArrayNew_ref = [.false._LKG, .true._LKG, .false._LKG]
#endif
            allocate(instance(0))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in ArrayTemplate sandwiched with two Replacements as new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            ArrayNew_ref = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is an empty array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is an empty array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is an empty array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is an empty array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is an empty array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is an empty array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is an empty array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is an empty array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is an empty array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = "A"
            pattern = ""
            Replacement = "XX"
            ArrayNew_ref = "XXAXX"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            ArrayTemplate = ["AA"]
            allocate(character(2,SKG) :: pattern(0))
            Replacement = ["XX"]
            ArrayNew_ref = ["XX", "AA", "XX"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            ArrayTemplate = [33_IKG]
            allocate(pattern(0))
            Replacement = [11_IKG]
            ArrayNew_ref = [11_IKG, 33_IKG, 11_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            ArrayTemplate = [33._CKG]
            allocate(pattern(0))
            Replacement = [11._CKG]
            ArrayNew_ref = [11._CKG, 33._CKG, 11._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            ArrayTemplate = [33._RKG]
            allocate(pattern(0))
            Replacement = [11._RKG]
            ArrayNew_ref = [11._RKG, 33._RKG, 11._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            ArrayTemplate = [.true._LKG]
            allocate(pattern(0))
            Replacement = [.false._LKG]
            ArrayNew_ref = [.false._LKG, .true._LKG, .false._LKG]
#endif
            instance = [1_IK, -1_IK]

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in ArrayTemplate sandwiched with two Replacements as new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty sorted unique array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty sorted unique array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty sorted unique array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty sorted unique array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty sorted unique array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty sorted unique array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty sorted unique array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty sorted unique array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty sorted unique array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = "A"
            pattern = ""
            Replacement = "XX"
            ArrayNew_ref = "XXAXX"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            ArrayTemplate = ["AA"]
            allocate(character(2,SKG) :: pattern(0))
            Replacement = ["XX"]
            ArrayNew_ref = ["XX", "AA", "XX"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            ArrayTemplate = [33_IKG]
            allocate(pattern(0))
            Replacement = [11_IKG]
            ArrayNew_ref = [11_IKG, 33_IKG, 11_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            ArrayTemplate = [33._CKG]
            allocate(pattern(0))
            Replacement = [11._CKG]
            ArrayNew_ref = [11._CKG, 33._CKG, 11._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            ArrayTemplate = [33._RKG]
            allocate(pattern(0))
            Replacement = [11._RKG]
            ArrayNew_ref = [11._RKG, 33._RKG, 11._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            ArrayTemplate = [.true._LK]
            allocate(pattern(0))
            Replacement = [.false._LK]
            ArrayNew_ref = [.false._LK, .true._LK, .false._LK]
#endif
            instance = [-1000_IK, 250_IK]

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in ArrayTemplate sandwiched with two Replacements as new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            ArrayNew_ref = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty insane array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty insane array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty insane array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty insane array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty insane array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty insane array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty insane array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty insane array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty insane array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = "A"
            pattern = ""
            Replacement = "XX"
            ArrayNew_ref = "AXX"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            ArrayTemplate = ["AA"]
            allocate(character(2,SKG) :: pattern(0))
            Replacement = ["XX"]
            ArrayNew_ref = ["AA", "XX"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            ArrayTemplate = [33_IKG]
            allocate(pattern(0))
            Replacement = [11_IKG]
            ArrayNew_ref = [33_IKG, 11_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            ArrayTemplate = [33._CKG]
            allocate(pattern(0))
            Replacement = [11._CKG]
            ArrayNew_ref = [33._CKG, 11._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            ArrayTemplate = [33._RKG]
            allocate(pattern(0))
            Replacement = [11._RKG]
            ArrayNew_ref = [33._RKG, 11._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            ArrayTemplate = [.true._LK]
            allocate(pattern(0))
            Replacement = [.false._LK]
            ArrayNew_ref = [.true._LK, .false._LK]
#endif

            Array = ArrayTemplate
            instance = [0_IK, 2_IK]
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in ArrayTemplate sandwiched with two Replacements as new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            instance = [0_IK, -1_IK]
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = ""
            pattern = "A"
            Replacement = "XX"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            allocate(character(2,SKG) :: ArrayTemplate(0))
            pattern = ["AA"]
            Replacement = ["XX"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            allocate(ArrayTemplate(0))
            pattern = [33_IKG]
            Replacement = [11_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            allocate(ArrayTemplate(0))
            pattern = [33._CKG]
            Replacement = [11._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            allocate(ArrayTemplate(0))
            pattern = [33._RKG]
            Replacement = [11._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            allocate(ArrayTemplate(0))
            pattern = [.true._LK]
            Replacement = [.false._LK]
#endif
            ArrayNew_ref = ArrayTemplate
            instance = [0_IK, -1_IK]

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in ArrayTemplate sandwiched with two Replacements as new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, with non-empty pattern and non-empty replacement, must result in the same array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = "AA"
            pattern = "A"
            Replacement = ""
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            ArrayTemplate = ["AA", "AA"]
            pattern = ["AA"]
            allocate(Replacement(0))
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            ArrayTemplate = [33_IKG, 33_IKG]
            pattern = [33_IKG]
            allocate(Replacement(0))
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            ArrayTemplate = [33._CKG, 33._CKG]
            pattern = [33._CKG]
            allocate(Replacement(0))
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            ArrayTemplate = [33._RKG, 33._RKG]
            pattern = [33._RKG]
            allocate(Replacement(0))
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            ArrayTemplate = [.true._LK]
            pattern = [.true._LK]
            allocate(Replacement(0))
#endif
            ArrayNew_ref = Replacement
            instance = [1_IK, -1_IK]

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and empty replacement, must result in an empty new array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = "AABBAACCAADDAAEE"
            pattern = "A"
            Replacement = "XX"
            ArrayNew_ref = "XXABBXXACCAADDXXXXEE"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            ArrayTemplate = ["AA", "AA", "BB", "BB", "AA", "AA", "CC", "CC", "AA", "AA", "DD", "DD", "AA", "AA", "EE", "EE"]
            pattern = ["AA"]
            Replacement = ["XX", "XX"]
            ArrayNew_ref = ["XX", "XX", "AA", "BB", "BB", "XX", "XX", "AA", "CC", "CC", "AA", "AA", "DD", "DD", "XX", "XX", "XX", "XX", "EE", "EE"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            ArrayTemplate = [1_IKG, 1_IKG, 2_IKG, 2_IKG, 1_IKG, 1_IKG, 3_IKG, 3_IKG, 1_IKG, 1_IKG, 4_IKG, 4_IKG, 1_IKG, 1_IKG, 5_IKG, 5_IKG]
            pattern = [1_IKG]
            Replacement = [0_IKG, 0_IKG]
            ArrayNew_ref = [0_IKG, 0_IKG, 1_IKG, 2_IKG, 2_IKG, 0_IKG, 0_IKG, 1_IKG, 3_IKG, 3_IKG, 1_IKG, 1_IKG, 4_IKG, 4_IKG, 0_IKG, 0_IKG, 0_IKG, 0_IKG, 5_IKG, 5_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            ArrayTemplate = [1._CKG, 1._CKG, 2._CKG, 2._CKG, 1._CKG, 1._CKG, 3._CKG, 3._CKG, 1._CKG, 1._CKG, 4._CKG, 4._CKG, 1._CKG, 1._CKG, 5._CKG, 5._CKG]
            pattern = [1._CKG]
            Replacement = [0._CKG, 0._CKG]
            ArrayNew_ref = [0._CKG, 0._CKG, 1._CKG, 2._CKG, 2._CKG, 0._CKG, 0._CKG, 1._CKG, 3._CKG, 3._CKG, 1._CKG, 1._CKG, 4._CKG, 4._CKG, 0._CKG, 0._CKG, 0._CKG, 0._CKG, 5._CKG, 5._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            ArrayTemplate = [1._RKG, 1._RKG, 2._RKG, 2._RKG, 1._RKG, 1._RKG, 3._RKG, 3._RKG, 1._RKG, 1._RKG, 4._RKG, 4._RKG, 1._RKG, 1._RKG, 5._RKG, 5._RKG]
            pattern = [1._RKG]
            Replacement = [0._RKG, 0._RKG]
            ArrayNew_ref = [0._RKG, 0._RKG, 1._RKG, 2._RKG, 2._RKG, 0._RKG, 0._RKG, 1._RKG, 3._RKG, 3._RKG, 1._RKG, 1._RKG, 4._RKG, 4._RKG, 0._RKG, 0._RKG, 0._RKG, 0._RKG, 5._RKG, 5._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            ArrayTemplate = [.true._LK, .true._LK, .false._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK]
            pattern = [.true._LK]
            Replacement = [.false._LK, .false._LK]
            ArrayNew_ref = [.false._LK, .false._LK, .true._LK, .false._LK, .false._LK, .false._LK, .false._LK, .true._LK, .false._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK, .false._LK, .false._LK, .false._LK, .false._LK, .false._LK, .false._LK]
#endif

            Array = ArrayTemplate
            instance = [1_IK, 3_IK, 7_IK, 8_IK]
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            instance = [1_IK, 3_IK, -2_IK, -1_IK]
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = "AABBAACCAADDAAEE"
            pattern = "A"
            Replacement = "XX"
            ArrayNew_ref = "XXXXBBXXXXCCXXXXDDXXXXEE"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            ArrayTemplate = ["AA", "AA", "BB", "BB", "AA", "AA", "CC", "CC", "AA", "AA", "DD", "DD", "AA", "AA", "EE", "EE"]
            pattern = ["AA"]
            Replacement = ["XX", "XX"]
            ArrayNew_ref = ["XX", "XX", "XX", "XX", "BB", "BB", "XX", "XX", "XX", "XX", "CC", "CC", "XX", "XX", "XX", "XX", "DD", "DD", "XX", "XX", "XX", "XX", "EE", "EE"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            ArrayTemplate = [1_IKG, 1_IKG, 2_IKG, 2_IKG, 1_IKG, 1_IKG, 3_IKG, 3_IKG, 1_IKG, 1_IKG, 4_IKG, 4_IKG, 1_IKG, 1_IKG, 5_IKG, 5_IKG]
            pattern = [1_IKG]
            Replacement = [0_IKG, 0_IKG]
            ArrayNew_ref = [0_IKG, 0_IKG, 0_IKG, 0_IKG, 2_IKG, 2_IKG, 0_IKG, 0_IKG, 0_IKG, 0_IKG, 3_IKG, 3_IKG, 0_IKG, 0_IKG, 0_IKG, 0_IKG, 4_IKG, 4_IKG, 0_IKG, 0_IKG, 0_IKG, 0_IKG, 5_IKG, 5_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            ArrayTemplate = [1._CKG, 1._CKG, 2._CKG, 2._CKG, 1._CKG, 1._CKG, 3._CKG, 3._CKG, 1._CKG, 1._CKG, 4._CKG, 4._CKG, 1._CKG, 1._CKG, 5._CKG, 5._CKG]
            pattern = [1._CKG]
            Replacement = [0._CKG, 0._CKG]
            ArrayNew_ref = [0._CKG, 0._CKG, 0._CKG, 0._CKG, 2._CKG, 2._CKG, 0._CKG, 0._CKG, 0._CKG, 0._CKG, 3._CKG, 3._CKG, 0._CKG, 0._CKG, 0._CKG, 0._CKG, 4._CKG, 4._CKG, 0._CKG, 0._CKG, 0._CKG, 0._CKG, 5._CKG, 5._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            ArrayTemplate = [1._RKG, 1._RKG, 2._RKG, 2._RKG, 1._RKG, 1._RKG, 3._RKG, 3._RKG, 1._RKG, 1._RKG, 4._RKG, 4._RKG, 1._RKG, 1._RKG, 5._RKG, 5._RKG]
            pattern = [1._RKG]
            Replacement = [0._RKG, 0._RKG]
            ArrayNew_ref = [0._RKG, 0._RKG, 0._RKG, 0._RKG, 2._RKG, 2._RKG, 0._RKG, 0._RKG, 0._RKG, 0._RKG, 3._RKG, 3._RKG, 0._RKG, 0._RKG, 0._RKG, 0._RKG, 4._RKG, 4._RKG, 0._RKG, 0._RKG, 0._RKG, 0._RKG, 5._RKG, 5._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            ArrayTemplate = [.true., .true., .false., .false., .true., .true., .false., .false., .true., .true., .false., .false., .true., .true., .false., .false.]
            pattern = [.true.]
            Replacement = [.false., .false.]
            ArrayNew_ref = [.false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false., .false.]
#endif
            instance = [-8_IK, -7_IK, -6_IK, -5_IK, -4_IK, -3_IK, -2_IK, -1_IK]

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            ! unsort and duplicate instance elements.
            call setShuffled(instance)
            instance = [instance, 1_IK]

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty non-unique unsorted mixed array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty non-unique unsorted mixed array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = "AABBAACCAADDAAEE"
            pattern = "Y"
            Replacement = "XX"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            ArrayTemplate = ["AA", "AA", "BB", "BB", "AA", "AA", "CC", "CC", "AA", "AA", "DD", "DD", "AA", "AA", "EE", "EE"]
            pattern = ["YY"]
            Replacement = ["XX", "XX"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            ArrayTemplate = [1_IKG, 1_IKG, 2_IKG, 2_IKG, 1_IKG, 1_IKG, 3_IKG, 3_IKG, 1_IKG, 1_IKG, 4_IKG, 4_IKG, 1_IKG, 1_IKG, 5_IKG, 5_IKG]
            pattern = [-1_IKG]
            Replacement = [0_IKG, 0_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            ArrayTemplate = [1._CKG, 1._CKG, 2._CKG, 2._CKG, 1._CKG, 1._CKG, 3._CKG, 3._CKG, 1._CKG, 1._CKG, 4._CKG, 4._CKG, 1._CKG, 1._CKG, 5._CKG, 5._CKG]
            pattern = [-1._CKG]
            Replacement = [0._CKG, 0._CKG]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            ArrayTemplate = [1._RKG, 1._RKG, 2._RKG, 2._RKG, 1._RKG, 1._RKG, 3._RKG, 3._RKG, 1._RKG, 1._RKG, 4._RKG, 4._RKG, 1._RKG, 1._RKG, 5._RKG, 5._RKG]
            pattern = [-1._RKG]
            Replacement = [0._RKG, 0._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            ArrayTemplate = [.true., .true., .true., .true., .true., .true., .true., .true., .true., .true., .true., .true., .true., .true., .true., .true.]
            pattern = [.false.]
            Replacement = [.false., .false.]
#endif
            instance = [-8_IK, -7_IK, -6_IK, -5_IK, -4_IK, -3_IK, -2_IK, -1_IK]

            ArrayNew_ref = ArrayTemplate
            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty mixed array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            ! unsort and duplicate instance elements.
            call setShuffled(instance)
            instance = [instance, 1_IK]

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty non-unique unsorted mixed array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty array, with non-empty matching pattern and non-empty replacement, must result in a properly-replaced new array when `instance` is a non-empty non-unique unsorted mixed array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getReplaced_D0_D0_D0_SK_ENABLED || setReplaced_D0_D0_D0_SK_ENABLED
            ArrayTemplate = "AA"
            Replacement = "XX"
#elif       getReplaced_D1_D1_D1_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            ArrayTemplate = ["AA"]
            Replacement = ["XX"]
#elif       getReplaced_D1_D1_D1_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            ArrayTemplate = [33_IKG]
            Replacement = [11_IKG, 11_IKG]
#elif       getReplaced_D1_D1_D1_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            ArrayTemplate = [(33._CKG,-33._CKG)]
            Replacement = [(11._CKG,-11._CKG), (11._CKG,-11._CKG)]
#elif       getReplaced_D1_D1_D1_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            ArrayTemplate = [33._RKG]
            Replacement = [11._RKG, 11._RKG]
#elif       getReplaced_D1_D1_D1_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            ArrayTemplate = [.true., .true.]
            Replacement = [.false., .false.]
#endif
            pattern = ArrayTemplate
            ArrayNew_ref = ArrayTemplate
            instance = [0_IK, 3_IK]

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a nonsensical array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a nonsensical array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a nonsensical array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a nonsensical array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a nonsensical array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a nonsensical array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a nonsensical array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a nonsensical array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a nonsensical array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            instance = [1_IK, 2_IK, -1_IK]
            ArrayNew_ref = Replacement

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a mixed but sensible array.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a mixed but sensible array with `sorted = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a mixed but sensible array with `sorted = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a mixed but sensible array with `unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a mixed but sensible array with `unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a mixed but sensible array with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a mixed but sensible array with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a mixed but sensible array with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            Array = ArrayTemplate
            call report(getReplacedEnabled, iseq, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An non-empty array that identical to the pattern with non-empty replacement, must result in the same array when `instance` is a mixed but sensible array with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else
            !%%%%%%%%%%%%%%%%%%%%%%%%
#error      "Unrecognized interface."
            !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(getReplacedEnabled, iseq, instance, sorted, unique)
            logical(LK) , intent(in)                        :: getReplacedEnabled
            logical(LK) , external  , optional              :: iseq
            integer(IK) , intent(in), optional, contiguous  :: instance(:)
            logical(LK) , intent(in), optional              :: sorted
            logical(LK) , intent(in), optional              :: unique
            !write(*,"(*(g0,:,', '))") "Array", Array
            !write(*,"(*(g0,:,', '))") "pattern", pattern
            !write(*,"(*(g0,:,', '))") "Replacement", Replacement
            !write(*,"(*(g0,:,', '))") "present(iseq)", present(iseq)
            !if (present(instance)) write(*,"(*(g0,:,', '))") "instance", instance
            !if (present(sorted)) write(*,"(*(g0,:,', '))") "sorted", sorted
            !if (present(unique)) write(*,"(*(g0,:,', '))") "unique", unique
            if (present(iseq) .and. present(instance)) then
#if             getReplaced_D1_D0_D0_ENABLED || getReplaced_D1_D0_D1_ENABLED || getReplaced_D1_D1_D0_ENABLED || \
                setReplaced_D1_D0_D0_ENABLED || setReplaced_D1_D0_D1_ENABLED || setReplaced_D1_D1_D0_ENABLED
                ArrayNew_ref = getReplaced(Array, [pattern], [Replacement], iseqVec, instance, sorted, unique)
#endif
                if (getReplacedEnabled) then
                    Array = getReplaced(Array, pattern, Replacement, iseq, instance, sorted, unique)
                else
                    call setReplaced(Array, pattern, Replacement, iseq, instance, sorted, unique)
                end if

            elseif (present(iseq) .and. .not. present(instance)) then
#if             getReplaced_D1_D0_D0_ENABLED || getReplaced_D1_D0_D1_ENABLED || getReplaced_D1_D1_D0_ENABLED || \
                setReplaced_D1_D0_D0_ENABLED || setReplaced_D1_D0_D1_ENABLED || setReplaced_D1_D1_D0_ENABLED
                ArrayNew_ref = getReplaced(Array, [pattern], [Replacement], iseqVec)
#endif
                if (getReplacedEnabled) then
                    Array = getReplaced(Array, pattern, Replacement, iseq)
                else
                    call setReplaced(Array, pattern, Replacement, iseq)
                end if
            elseif (.not. present(iseq) .and. present(instance)) then
#if             getReplaced_D1_D0_D0_ENABLED || getReplaced_D1_D0_D1_ENABLED || getReplaced_D1_D1_D0_ENABLED || \
                setReplaced_D1_D0_D0_ENABLED || setReplaced_D1_D0_D1_ENABLED || setReplaced_D1_D1_D0_ENABLED
                ArrayNew_ref = getReplaced(Array, [pattern], [Replacement], instance, sorted, unique)
#endif
                if (getReplacedEnabled) then
                    Array = getReplaced(Array, pattern, Replacement, instance, sorted, unique)
                else
                    call setReplaced(Array, pattern, Replacement, instance, sorted, unique)
                end if
            else
#if             getReplaced_D1_D0_D0_ENABLED || getReplaced_D1_D0_D1_ENABLED || getReplaced_D1_D1_D0_ENABLED || \
                setReplaced_D1_D0_D0_ENABLED || setReplaced_D1_D0_D1_ENABLED || setReplaced_D1_D1_D0_ENABLED
                ArrayNew_ref = getReplaced(Array, [pattern], [Replacement])
#endif
                if (getReplacedEnabled) then
                    Array = getReplaced(Array, pattern, Replacement)
                else
                    call setReplaced(Array, pattern, Replacement)
                end if
            end if

            assertion = assertion .and. GET_SIZE(Array, kind = IK) == GET_SIZE(ArrayNew_ref, kind = IK)
            call setOutput(getReplacedEnabled, iseq, instance, sorted, unique)
            call test%assert(assertion, PROCEDURE_NAME//SK_": The length of the output array must match the length of the reference array.", int(__LINE__, IK))

            assertion = assertion .and. ALL(Array IS_EQUAL ArrayNew_ref)
            call setOutput(getReplacedEnabled, iseq, instance, sorted, unique)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine setOutput(getReplacedEnabled, iseq, instance, sorted, unique)
            logical(LK) , intent(in)                        :: getReplacedEnabled
            logical(LK) , external  , optional              :: iseq
            integer(IK) , intent(in), optional, contiguous  :: instance(:)
            logical(LK) , intent(in), optional              :: sorted
            logical(LK) , intent(in), optional              :: unique
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Array              ", Array
                write(test%disp%unit,"(*(g0,:,', '))") "ArrayNew_ref       ", ArrayNew_ref
                write(test%disp%unit,"(*(g0,:,', '))") "ArrayTemplate      ", ArrayTemplate
                write(test%disp%unit,"(*(g0,:,', '))") "size(Array)        ", GET_SIZE(Array, kind = IK)
                write(test%disp%unit,"(*(g0,:,', '))") "size(ArrayNew_ref) ", GET_SIZE(ArrayNew_ref, kind = IK)
                write(test%disp%unit,"(*(g0,:,', '))") "size(ArrayTemplate)", GET_SIZE(ArrayTemplate, kind = IK)
                write(test%disp%unit,"(*(g0,:,', '))") "Replacement        ", Replacement
                write(test%disp%unit,"(*(g0,:,', '))") "pattern            ", pattern
                write(test%disp%unit,"(*(g0,:,', '))") "present(instance)  ", present(instance)
                write(test%disp%unit,"(*(g0,:,', '))") "present(unique)    ", present(unique)
                write(test%disp%unit,"(*(g0,:,', '))") "present(sorted)    ", present(sorted)
                write(test%disp%unit,"(*(g0,:,', '))") "present(iseq)      ", present(iseq)
                if (present(instance)) then
                write(test%disp%unit,"(*(g0,:,', '))") "instance           ", instance
                end if
                if (present(sorted)) then
                write(test%disp%unit,"(*(g0,:,', '))") "sorted             ", sorted
                end if
                if (present(unique)) then
                write(test%disp%unit,"(*(g0,:,', '))") "unique             ", unique
                end if
                if (present(instance)) then
                write(test%disp%unit,"(*(g0,:,', '))") "getReplacedEnabled ", getReplacedEnabled
                end if
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            implicit none
            if (allocated(ArrayTemplate)) deallocate(ArrayTemplate)
#if         getReplaced_D1_D1_D0_ENABLED || setReplaced_D1_D1_D0_ENABLED || \
            getReplaced_D1_D1_D1_ENABLED || setReplaced_D1_D1_D1_ENABLED
            if (allocated(pattern)) deallocate(pattern)
#endif
#if         getReplaced_D1_D0_D1_ENABLED || setReplaced_D1_D0_D1_ENABLED || \
            getReplaced_D1_D1_D1_ENABLED || setReplaced_D1_D1_D1_ENABLED
            if (allocated(Replacement)) deallocate(Replacement)
#endif
            if (allocated(ArrayNew_ref)) deallocate(ArrayNew_ref)
            if (allocated(instance)) deallocate(instance)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getReplaced_D0_D0_D0_ENABLED || setReplaced_D0_D0_D0_ENABLED
        function iseq(Segment, pattern) result(equivalent)
            character(*,SKG), intent(in) :: pattern, Segment
            logical(LK) :: equivalent
            equivalent = Segment == pattern
        end function
#elif   getReplaced_D1_D0_D0_ENABLED || getReplaced_D1_D0_D1_ENABLED || \
        setReplaced_D1_D0_D0_ENABLED || setReplaced_D1_D0_D1_ENABLED
        function iseq(Segment, pattern) result(equivalent)
#if         getReplaced_D1_D0_D0_SK_ENABLED || getReplaced_D1_D0_D1_SK_ENABLED || \
            setReplaced_D1_D0_D0_SK_ENABLED || setReplaced_D1_D0_D1_SK_ENABLED
            character(*,SKG), intent(in)    :: Segment, pattern
#elif       getReplaced_D1_D0_D0_IK_ENABLED || getReplaced_D1_D0_D1_IK_ENABLED || \
            setReplaced_D1_D0_D0_IK_ENABLED || setReplaced_D1_D0_D1_IK_ENABLED
            integer(IKG)    , intent(in)    :: Segment, pattern
#elif       getReplaced_D1_D0_D0_LK_ENABLED || getReplaced_D1_D0_D1_LK_ENABLED || \
            setReplaced_D1_D0_D0_LK_ENABLED || setReplaced_D1_D0_D1_LK_ENABLED
            logical(LKG)    , intent(in)    :: Segment, pattern
#elif       getReplaced_D1_D0_D0_CK_ENABLED || getReplaced_D1_D0_D1_CK_ENABLED || \
            setReplaced_D1_D0_D0_CK_ENABLED || setReplaced_D1_D0_D1_CK_ENABLED
            complex(CKG)    , intent(in)    :: Segment, pattern
#elif       getReplaced_D1_D0_D0_RK_ENABLED || getReplaced_D1_D0_D1_RK_ENABLED || \
            setReplaced_D1_D0_D0_RK_ENABLED || setReplaced_D1_D0_D1_RK_ENABLED
            real(RKG)       , intent(in)    :: Segment, pattern
#endif
            logical(LK) :: equivalent
            equivalent = Segment IS_EQUAL pattern
        end function
#elif   getReplaced_D1_D1_D0_ENABLED || getReplaced_D1_D1_D1_ENABLED || \
        setReplaced_D1_D1_D0_ENABLED || setReplaced_D1_D1_D1_ENABLED
        function iseq(Segment, pattern, lenPattern) result(equivalent)
            logical(LK)             :: equivalent
            integer(IK), intent(in) :: lenPattern
#if         getReplaced_D1_D1_D0_SK_ENABLED || getReplaced_D1_D1_D1_SK_ENABLED || \
            setReplaced_D1_D1_D0_SK_ENABLED || setReplaced_D1_D1_D1_SK_ENABLED
            character(*,SKG), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       getReplaced_D1_D1_D0_IK_ENABLED || getReplaced_D1_D1_D1_IK_ENABLED || \
            setReplaced_D1_D1_D0_IK_ENABLED || setReplaced_D1_D1_D1_IK_ENABLED
            integer(IKG)    , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       getReplaced_D1_D1_D0_LK_ENABLED || getReplaced_D1_D1_D1_LK_ENABLED || \
            setReplaced_D1_D1_D0_LK_ENABLED || setReplaced_D1_D1_D1_LK_ENABLED
            logical(LKG)    , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       getReplaced_D1_D1_D0_CK_ENABLED || getReplaced_D1_D1_D1_CK_ENABLED || \
            setReplaced_D1_D1_D0_CK_ENABLED || setReplaced_D1_D1_D1_CK_ENABLED
            complex(CKG)    , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       getReplaced_D1_D1_D0_RK_ENABLED || getReplaced_D1_D1_D1_RK_ENABLED || \
            setReplaced_D1_D1_D0_RK_ENABLED || setReplaced_D1_D1_D1_RK_ENABLED
            real(RKG)       , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#endif
            equivalent = all(Segment IS_EQUAL pattern)
        end function
#else
#error  "Unrecognized interface."
#endif

#if     getReplaced_D1_D0_D0_ENABLED || getReplaced_D1_D0_D1_ENABLED || getReplaced_D1_D1_D0_ENABLED || \
        setReplaced_D1_D0_D0_ENABLED || setReplaced_D1_D0_D1_ENABLED || setReplaced_D1_D1_D0_ENABLED
        function iseqVec(Segment, pattern, lenPattern) result(equivalent)
            logical(LK)             :: equivalent
            integer(IK), intent(in) :: lenPattern
#if         getReplaced_D1_D0_D0_SK_ENABLED || getReplaced_D1_D0_D1_SK_ENABLED || getReplaced_D1_D1_D0_SK_ENABLED || \
            setReplaced_D1_D0_D0_SK_ENABLED || setReplaced_D1_D0_D1_SK_ENABLED || setReplaced_D1_D1_D0_SK_ENABLED
            character(*,SKG), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       getReplaced_D1_D0_D0_IK_ENABLED || getReplaced_D1_D0_D1_IK_ENABLED || getReplaced_D1_D1_D0_IK_ENABLED || \
            setReplaced_D1_D0_D0_IK_ENABLED || setReplaced_D1_D0_D1_IK_ENABLED || setReplaced_D1_D1_D0_IK_ENABLED
            integer(IKG)    , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       getReplaced_D1_D0_D0_LK_ENABLED || getReplaced_D1_D0_D1_LK_ENABLED || getReplaced_D1_D1_D0_LK_ENABLED || \
            setReplaced_D1_D0_D0_LK_ENABLED || setReplaced_D1_D0_D1_LK_ENABLED || setReplaced_D1_D1_D0_LK_ENABLED
            logical(LKG)    , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       getReplaced_D1_D0_D0_CK_ENABLED || getReplaced_D1_D0_D1_CK_ENABLED || getReplaced_D1_D1_D0_CK_ENABLED || \
            setReplaced_D1_D0_D0_CK_ENABLED || setReplaced_D1_D0_D1_CK_ENABLED || setReplaced_D1_D1_D0_CK_ENABLED
            complex(CKG)    , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       getReplaced_D1_D0_D0_RK_ENABLED || getReplaced_D1_D0_D1_RK_ENABLED || getReplaced_D1_D1_D0_RK_ENABLED || \
            setReplaced_D1_D0_D0_RK_ENABLED || setReplaced_D1_D0_D1_RK_ENABLED || setReplaced_D1_D1_D0_RK_ENABLED
            real(RKG)       , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#endif
            equivalent = all(Segment IS_EQUAL pattern)
        end function
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  IS_EQUAL
#undef  GET_SIZE
#undef  ALL
