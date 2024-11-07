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
!>  [getLoc](@ref pm_arrayFind::getLoc).
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
        character(*,SK), parameter :: PROCEDURE_NAME = "@setLoc/getLoc()"

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     loc_ENABLED && D1_D0_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D1_D0_ENABLED
        character(2,SKG), dimension(:), allocatable :: array
        character(2,SKG)                             :: pattern
#elif   IK_ENABLED && D1_D0_ENABLED
        integer(IKG)    , dimension(:), allocatable :: array
        integer(IKG)                                :: pattern
#elif   LK_ENABLED && D1_D0_ENABLED
        logical(LKG)    , dimension(:), allocatable :: array
        logical(LKG)                                :: pattern
#elif   CK_ENABLED && D1_D0_ENABLED
        complex(CKG)    , dimension(:), allocatable :: array
        complex(CKG)                                :: pattern
#elif   RK_ENABLED && D1_D0_ENABLED
        real(RKG)       , dimension(:), allocatable :: array
        real(RKG)                                   :: pattern
#else
#error  "Unrecognized interface."
#endif
        integer(IK) , allocatable   :: loc(:)
        integer(IK) , allocatable   :: instance(:)
        integer(IK) , allocatable   :: loc_ref(:)
        logical(LK)                 :: setLocEnabled

        assertion = .true._LK

        setLocEnabled = .false._LK
        call runTestsWith()
        call runTestsWith(iseq = iseq)

        setLocEnabled = .true._LK
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
            if (allocated(loc_ref)) deallocate(loc_ref)
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
            loc_ref = [integer(IK) ::]
            allocate(instance(0))
            allocate(array(0))

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `blindness = 1_IK, sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `blindness = 1_IK, sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `blindness = 1_IK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `blindness = 1_IK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `blindness = 1_IK, sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `blindness = 1_IK, sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `blindness = 1_IK, sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting loc with `instance` with `blindness = 1_IK, sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            loc_ref= [1_IK, 2_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `blindness = 1_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield an `loc` of elements of `array` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            pattern = "BB"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = 0_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (0._CKG,-0._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = 0._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            pattern = .true._LK
#endif
            instance = [1_IK, 2_IK]
            allocate(loc_ref(0))

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` with `blindness = 1_IK` setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `blindness = 1_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are not the scalar `pattern` must yield an empty `loc` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            loc_ref = [1_IK, 2_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `blindness = 1_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `blindness = 1_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `blindness = 1_IK`, `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield the full index of `array` when `instance` is non-empty and out-of-bound with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            loc_ref = [1_IK, 2_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `blindness = 1_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `blindness = 1_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `blindness = 1_IK`, `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a full `loc` when `instance` is non-empty and mixed with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            loc_ref = [1_IK, 1_IK, 2_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            loc_ref = [1_IK, 2_IK, 1_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all scalar `pattern` must yield a redundant full `loc` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "BB", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 2_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 2._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .true._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [1_IK, -2_IK]
            loc_ref = [1_IK, 1_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc = [1,1]` when `instance = [1_IK, -2_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "BB", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 2_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 2._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .true._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [-2_IK]
            loc_ref = [1_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive with `blindness = 2_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `blindness = 2_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `blindness = 2_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `blindness = 2_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `blindness = 2_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the beginning and the end matching the scalar `pattern` must yield an `loc_ref = [1_IK]`  when `instance = [-2_IK]` is non-empty but incomprehensive and positive, with `blindness = 2_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["BB", "AA", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [2_IKG, 1_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(2._CKG,-2._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [2._RKG, 1._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.true._LK, .false._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [-1_IK, 3_IK]
            loc_ref = [3_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [2_IK]

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            instance = [-1_IK]

            loc_ref = [2_IK]
            call report(__LINE__, iseq, instance = instance, blindness = 5_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed with `blindness = 5_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 5_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 5_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 5_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 5_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 5_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 5_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 5_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 5_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 5_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 5_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 5_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 5_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 5_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 3 with the last two elements matching the scalar `pattern` must yield `loc = [3_IK]` when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 5_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            pattern = array(1)
            instance = [-1_IK, 3_IK]
            loc_ref = [1_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed with `blindness = 4_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 4_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 4_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 4_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 4_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 4_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 4_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield loc = [1_IK] when `instance = [-1_IK, 3_IK]` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 4_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            pattern = array(1)
            instance = [10_IK, 3_IK]
            loc_ref = [integer(IK)::]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [3_IK, 10_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 7_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 7_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 7_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 7_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 7_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 7_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 7_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 7_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 7_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 7_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 7_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 7_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 7_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 7_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            allocate(character(2,SKG) :: array(0))
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
            instance = [10_IK, 3_IK]
            loc_ref = [integer(IK)::]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [3_IK, 10_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 10_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 10_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 10_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 10_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 10_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 10_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 10_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield an empty `loc` that is the same as the `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 10_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "BB", "CC", "AA"]
            pattern = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 2_IKG, 3_IKG, 1_IKG]
            pattern = 1_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG)]
            pattern = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 2._RKG, 3._RKG, 1._RKG]
            pattern = 1._RKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .true._LK, .true._LK, .false._LK]
            pattern = .false._LK
#endif
            instance = [1_IK, -1_IK]
            loc_ref = [1_IK, 4_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 3_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 3_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 3_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 3_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 3_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 3_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 3_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 3_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 3_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 3_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 3_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 3_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 3_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 3_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield an `loc = [1_IK, 4_IK]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 3_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            instance = [10_IK, 3_IK]
            loc_ref = [integer(IK)::]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [3_IK, 10_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 1_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 1_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield an empty `loc` of size one whose value is `array` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            loc_ref = [2_IK, 4_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` with `blindness = 2_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `blindness = 2_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `blindness = 2_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `blindness = 2_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `blindness = 2_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `blindness = 2_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `blindness = 2_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc = [1_IK, 5_IK]` when `instance` is ordered, positive, non-empty but sensible, with `blindness = 2_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

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
            instance = [10_IK, 3_IK]
            loc_ref = [integer(IK)::]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [3_IK, 10_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield an empty `loc` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        subroutine report   ( line & ! LCOV_EXCL_LINE
                            , iseq & ! LCOV_EXCL_LINE
                            , instance & ! LCOV_EXCL_LINE
                            , blindness & ! LCOV_EXCL_LINE
                            , positive & ! LCOV_EXCL_LINE
                            , sorted & ! LCOV_EXCL_LINE
                            )
            integer, intent(in) :: line
            logical(LK) , external  , optional              :: iseq
            integer(IK) , intent(in), optional, contiguous  :: instance(:)
            integer(IK) , intent(in), optional              :: blindness
            logical(LK) , intent(in), optional              :: positive
            logical(LK) , intent(in), optional              :: sorted
            logical(LK) :: sorted_def, positive_def
            integer(IK) :: blindness_def, nloc

            sorted_def = .false._LK; if (present(sorted)) sorted_def = sorted
            positive_def = .false._LK; if (present(positive)) positive_def = positive
            blindness_def = 1_IK; if (present(blindness)) blindness_def = blindness

            if (setLocEnabled) then
                if (present(iseq) .and. present(instance)) then
                    call setLoc(loc, nloc, array, pattern, iseq, instance, sorted_def, positive_def, blindness_def)
                elseif (present(instance)) then
                    call setLoc(loc, nloc, array, pattern, instance, sorted_def, positive_def, blindness_def)
                elseif (present(iseq)) then
                    call setLoc(loc, nloc, array, pattern, iseq, blindness_def)
                else
                    call setLoc(loc, nloc, array, pattern, blindness_def)
                end if
            else
                if (present(iseq) .and. present(instance) .and. present(sorted) .and. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, instance = instance, sorted = sorted, positive = positive, blindness = blindness)
                elseif (present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, instance = instance, sorted = sorted, blindness = blindness)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, instance = instance, positive = positive, blindness = blindness)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, instance = instance, blindness = blindness)
                elseif (present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, blindness = blindness)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. present(positive)) then
                    loc = getLoc(array, pattern, instance = instance, sorted = sorted, positive = positive, blindness = blindness)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, instance = instance, sorted = sorted, blindness = blindness)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(positive)) then
                    loc = getLoc(array, pattern, instance = instance, positive = positive, blindness = blindness)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, instance = instance, blindness = blindness)
                elseif (.not. present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, blindness = blindness)
                else
                    error stop PROCEDURE_NAME//SK_": Unrecognized interface in testing." ! LCOV_EXCL_LINE
                end if
                nloc = size(loc, 1, IK)
            end if

            ! Report test results if needed.

            !write(*,*) setLocEnabled, present(instance), present(sorted), present(positive)
            !write(*,*) array
            !write(*,*) pattern
            !write(*,*) loc
            !write(*,*) loc_ref

            assertion = assertion .and. nloc == size(loc_ref, 1, IK)
            if (assertion .and. 0_IK < nloc) assertion = assertion .and. all(loc(1 : nloc) == loc_ref)

            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("array")
                call test%disp%show( array )
                call test%disp%show("pattern")
                call test%disp%show( pattern )
                call test%disp%show("nloc")
                call test%disp%show( nloc )
                call test%disp%show("loc_ref")
                call test%disp%show( loc_ref )
                call test%disp%show("loc(1:nloc)")
                call test%disp%show( loc(1:nloc) )
                call test%disp%show("present(iseq)")
                call test%disp%show( present(iseq) )
                call test%disp%show("present(instance)")
                call test%disp%show( present(instance) )
                if (present(instance)) then
                    call test%disp%show("instance")
                    call test%disp%show( instance )
                end if
                call test%disp%show("present(positive)")
                call test%disp%show( present(positive) )
                if (present(positive)) then
                    call test%disp%show("positive")
                    call test%disp%show( positive )
                end if
                call test%disp%show("present(sorted)")
                call test%disp%show( present(sorted) )
                if (present(sorted)) then
                    call test%disp%show("sorted")
                    call test%disp%show( sorted )
                end if
                call test%disp%show("present(blindness)")
                call test%disp%show( present(blindness) )
                if (present(blindness)) then
                    call test%disp%show("blindness")
                    call test%disp%show( blindness )
                end if
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if

            if (nloc /= size(loc_ref, 1, IK)) call test%assert(assertion, PROCEDURE_NAME//SK_": The condition `nloc = size(loc_ref) must hold.", int(line, IK))

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   loc_ENABLED && (D0_D0_ENABLED || D1_D1_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_D0_ENABLED
#define ALL
        character(:,SKG), allocatable :: array, pattern
#elif   SK_ENABLED && D1_D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: array, pattern
#elif   IK_ENABLED && D1_D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: array, pattern
#elif   LK_ENABLED && D1_D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: array, pattern
#elif   CK_ENABLED && D1_D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: array, pattern
#elif   RK_ENABLED && D1_D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: array, pattern
#else
#error "Unrecognized interface."
#endif
        integer(IK) , allocatable   :: loc(:)
        integer(IK) , allocatable   :: instance(:)
        integer(IK) , allocatable   :: loc_ref(:)
        logical(LK)                 :: setLocEnabled

        assertion = .true._LK

        setLocEnabled = .false._LK
        call runTestsWith()
        call runTestsWith(iseq = iseq)

        setLocEnabled = .true._LK
        call runTestsWith()
        call runTestsWith(iseq = iseq)

    contains

        function iseq(Segment, pattern & ! LCOV_EXCL_LINE
#if     D1_D1_ENABLED
        , lenPattern & ! LCOV_EXCL_LINE
#endif
        ) result(equivalent)
#if         D1_D1_ENABLED
            integer(IK)     , intent(in)    :: lenPattern
#endif
#if         SK_ENABLED && D0_D0_ENABLED
            character(*, SK), intent(in) :: segment, pattern
#elif       SK_ENABLED && D1_D1_ENABLED
            character(*, SK), intent(in) :: Segment(lenPattern), pattern(lenPattern)
#elif       IK_ENABLED && D1_D1_ENABLED
            integer(IKG)    , intent(in) :: Segment(lenPattern), pattern(lenPattern)
#elif       CK_ENABLED && D1_D1_ENABLED
            complex(CKG)    , intent(in) :: Segment(lenPattern), pattern(lenPattern)
#elif       RK_ENABLED && D1_D1_ENABLED
            real(RKG)       , intent(in) :: Segment(lenPattern), pattern(lenPattern)
#elif       LK_ENABLED && D1_D1_ENABLED
            logical(LKG)    , intent(in) :: Segment(lenPattern), pattern(lenPattern)
#endif
            logical(LK) :: equivalent
            equivalent = ALL(pattern IS_EQUAL Segment)
        end function

        subroutine reset()
            if (allocated(array)) deallocate(array)
            if (allocated(pattern)) deallocate(pattern)
            if (allocated(instance)) deallocate(instance)
            if (allocated(loc_ref)) deallocate(loc_ref)
        end subroutine reset

        subroutine runTestsWith(iseq)
            logical(LK), external, optional  :: iseq

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            pattern = " "
            allocate(character(0,SKG) :: array)
#elif       SK_ENABLED && D1_D1_ENABLED
            pattern = [" "]
            allocate(character(2,SKG) :: array(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            pattern = [1_IKG]
            allocate(array(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            pattern = [1._CKG]
            allocate(array(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            pattern = [1._RKG]
            allocate(array(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            pattern = [.false._LK]
            allocate(array(0))
#endif
            allocate(instance(0))
            allocate(loc_ref(0))

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `blindness = 20_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `blindness = 20_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `blindness = 20_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `blindness = 20_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `blindness = 20_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `blindness = 20_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `blindness = 20_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `blindness = 20_IK`, `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `blindness = 20_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty `array` has empty resulting `loc` with `instance` with `blindness = 20_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            pattern = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            pattern = [.false._LK]
#endif
            instance = [1_IK, 2_IK]
            loc_ref = [1_IK, 2_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `blindness = 1_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            pattern = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            pattern = [.false._LK]
#endif
            instance = [1_IK, 2_IK, 3_IK]
            loc_ref = [1_IK, 2_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `blindness = 1_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `blindness = 1_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty and redundant with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            pattern = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            pattern = [.false._LK]
#endif
            instance = [1_IK, -1_IK]
            loc_ref = [1_IK, 2_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield loc = [1_IK, 2_IK] when `instance` is non-empty and mixed with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            pattern = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            pattern = [.false._LK]
#endif
            instance = [1_IK, -2_IK, -1_IK]
            loc_ref = [1_IK, 1_IK, 2_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the single-element vector `pattern` must yield `loc = [1_IK, 1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA", "AA"]
            pattern = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG, 1_IKG]
            pattern = [1_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG, 1._RKG]
            pattern = [1._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK, .false._LK]
            pattern = [.false._LK, .false._LK]
#endif
            instance = [1_IK, 2_IK]
            loc_ref = [1_IK, 2_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK]` when `instance` is non-empty, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA", "AA"]
            pattern = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG, 1_IKG]
            pattern = [1_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG, 1._RKG]
            pattern = [1._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK, .false._LK]
            pattern = [.false._LK, .false._LK]
#endif
            instance = [1_IK, -1_IK, -2_IK, 0_IK]
            loc_ref = [1_IK, 2_IK, 1_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK, 1_IK]` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK, 1_IK]` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK, 1_IK]` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK, 1_IK]` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK, 1_IK]` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK, 1_IK]` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK, 1_IK]` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` whose elements are all the double-element vector `pattern` must yield `loc = [1_IK, 2_IK, 1_IK]` when `instance` is non-empty, redundant, non-positive, unsorted, and mixed, with `blindness = 1_IK`, `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AABBAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA", "BB", "BB", "AA", "AA"]
            pattern = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG, 2_IKG, 2_IKG, 1_IKG, 1_IKG]
            pattern = [1_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG, 2._RKG, 2._RKG, 1._RKG, 1._RKG]
            pattern = [1._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK]
            pattern = [.false._LK, .false._LK]
#endif
            instance = [1_IK, -2_IK, -1_IK]
            loc_ref = [1_IK, 1_IK, 5_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK, 1_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AABBAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA", "BB", "BB", "AA", "AA"]
            pattern = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG, 2_IKG, 2_IKG, 1_IKG, 1_IKG]
            pattern = [1_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG, 2._RKG, 2._RKG, 1._RKG, 1._RKG]
            pattern = [1._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK]
            pattern = [.false._LK, .false._LK]
#endif
            instance = [-1_IK]
            loc_ref = [5_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [5_IK]` when `instance` is non-empty but incomprehensive, non-positive with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [5_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [5_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [5_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [5_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [5_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            loc_ref = [1_IK]

            call report(__LINE__, iseq, instance = instance, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK]` when `instance` is non-empty but incomprehensive, non-positive with `blindness = 6_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `blindness = 6_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `blindness = 6_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `blindness = 6_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `blindness = 6_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An `array` of length 5 with the beginning and the end matching the double-element vector `pattern` must yield an `loc = [1_IK]` when `instance` is non-empty but incomprehensive, non-positive, with `blindness = 6_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK]
#endif
            pattern = array
            instance = [-1_IK, 3_IK]
            loc_ref = [1_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 1_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "BBAAAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["BB", "BB", "AA", "AA", "AA", "AA"]
            pattern = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [2_IKG, 2_IKG, 1_IKG, 1_IKG, 1_IKG, 1_IKG]
            pattern = [1_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(2._CKG,-2._CKG), (2._CKG,-2._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [2._RKG, 2._RKG, 1._RKG, 1._RKG, 1._RKG, 1._RKG]
            pattern = [1._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.true._LK, .true._LK, .false._LK, .false._LK, .false._LK, .false._LK]
            pattern = [.false._LK, .false._LK]
#endif
            instance = [-1_IK, 3_IK]
            loc_ref = [5_IK, 5_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [3, 3]

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            instance = [-1_IK, 3_IK]
            loc_ref = [5_IK]

            call report(__LINE__, iseq, instance = instance, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed with `blindness = 2_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [2]

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with the last two elements matching the single-element vector `pattern` must yield `loc = [5_IK, 5_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK]
#endif
            pattern = array
            instance = [4_IK, 3_IK]
            loc_ref  = [integer(IK)::]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 1_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc  = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with `blindness = 1_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            allocate(character(0,SKG) :: array)
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            allocate(character(2,SKG) :: array(0))
            pattern = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            allocate(array(0))
            pattern = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            allocate(array(0))
            pattern = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            allocate(array(0))
            pattern = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            allocate(array(0))
            pattern = [.false._LK]
#endif
            loc_ref = [integer(IK)::]
            instance = [30_IK, 3_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 4_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 4_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 4_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 4_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 4_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 4_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 4_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 4_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 4_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 4_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 4_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 4_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 4_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 4_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `blindness = 4_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AABBCCAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA", "BB", "BB", "CC", "CC", "AA", "AA"]
            pattern = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG, 2_IKG, 2_IKG, 3_IKG, 3_IKG, 1_IKG, 1_IKG]
            pattern = [1_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG, 2._RKG, 2._RKG, 3._RKG, 3._RKG, 1._RKG, 1._RKG]
            pattern = [1._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK, .true._LK, .true._LK, .true._LK, .true._LK, .false._LK, .false._LK]
            pattern = [.false._LK, .false._LK]
#endif
            instance = [1_IK, -1_IK]
            loc_ref = [1_IK, 7_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, 2_IK]` with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, 2_IK]` with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, -1_IK]` with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, 2_IK]` with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 6_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` occurrence at the beginning and the end of the `array` must yield `loc = [1_IK, 7_IK]` when `instance = [1_IK, 2_IK]` with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AAAA"
            pattern = "XX"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            pattern = ["XX"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = [0_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(0._CKG,0._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = [0._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            pattern = [.true._LK]
#endif
            instance = [10_IK, 3_IK]
            loc_ref = [integer(IK)::]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [3_IK, 10_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` with `blindness = 2_IK`, setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, `blindness = 2_IK`, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, `blindness = 2_IK`, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, `blindness = 2_IK`, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, `blindness = 2_IK`, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, `blindness = 2_IK`, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, `blindness = 2_IK`, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, `blindness = 2_IK`, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, `blindness = 2_IK`, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, `blindness = 2_IK`, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, `blindness = 2_IK`, with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, `blindness = 2_IK`, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, `blindness = 2_IK`, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, `blindness = 2_IK`, with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `pattern` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty, and sensible, `blindness = 2_IK`, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AXAXA"
            pattern = "X"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "XX", "AA", "XX", "AA"]
            pattern = ["XX"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 0_IKG, 1_IKG, 0_IKG, 1_IKG]
            pattern = [0_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG)]
            pattern = [(0._CKG,0._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 0._RKG, 1._RKG, 0._RKG, 1._RKG]
            pattern = [0._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK]
            pattern = [.true._LK]
#endif
            loc_ref = [2_IK, 4_IK]

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance = [1_IK, 2_IK]` with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance = [1_IK, 2_IK]` with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            instance = [1_IK, -1_IK]
            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [1_IK, 2_IK]
            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance = [1_IK, 2_IK]` with `blindness = 2_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `pattern` occurrences must yield `loc_ref = [2_IK, 4_IK]` when `instance = [1_IK, 2_IK]` with `blindness = 2_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AAA"
            pattern = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            pattern = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = [1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = [1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            pattern = [.false._LK]
#endif
            instance = [10_IK, 3_IK]
            loc_ref = [integer(IK)::]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            instance = [3_IK, 10_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `blindness = 2_IK`, `sorted = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `blindness = 2_IK`, `sorted = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `blindness = 2_IK`, `positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `blindness = 2_IK`, `positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `blindness = 2_IK`, `sorted = .true._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `blindness = 2_IK`, `sorted = .true._LK, positive = .true.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `pattern` must yield `loc = [integer(IK)::]` when `instance = [0_IK, 3_IK]` with `blindness = 2_IK`, `sorted = .false._LK, positive = .false.` with setLocEnabled = "//getStr(setLocEnabled)//".", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            pattern = ""
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            allocate(character(2,SKG) :: pattern(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            allocate(pattern(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [1._CKG, 1._CKG]
            allocate(pattern(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            allocate(pattern(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            allocate(pattern(0))
#endif
            loc_ref = [integer(IK)::]
            allocate(instance(0))

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is empty.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is empty with `sorted = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is empty with `sorted = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is empty with `positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is empty with `positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is empty with `sorted = .true._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is empty with `sorted = .false._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is empty with `sorted = .true._LK, positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is empty with `sorted = .false._LK, positive = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            pattern = ""
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            allocate(character(2,SKG) :: pattern(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            allocate(pattern(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [1._CKG, 1._CKG]
            allocate(pattern(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            allocate(pattern(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            allocate(pattern(0))
#endif
            loc_ref = [integer(IK)::]
            instance = [1_IK]

            call report(__LINE__, iseq)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `sorted = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `sorted = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `sorted = .true._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `sorted = .false._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `sorted = .true._LK, positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `sorted = .false._LK, positive = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` with `blindness  = 2_IK`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `blindness  = 2_IK`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `blindness  = 2_IK`, `sorted = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `blindness  = 2_IK`, `sorted = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `blindness  = 2_IK`, `positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `blindness  = 2_IK`, `positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `blindness  = 2_IK`, `sorted = .true._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .true._LK, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `blindness  = 2_IK`, `sorted = .false._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `blindness  = 2_IK`, `sorted = .true._LK, positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness  = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `pattern` must yield `loc = [integer(IK)::]` when `instance` is non-empty with `blindness  = 2_IK`, `sorted = .false._LK, positive = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
#endif
            pattern = array
            instance = [-1_IK, 3_IK]
            loc_ref = [1_IK]

            call report(__LINE__, iseq, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .true._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `sorted = .false._LK, positive = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed with `blindness = 2_IK`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, with `sorted = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, with `sorted = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, with `positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, with `positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, with `sorted = .true._LK, positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, with `sorted = .true._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [1_IK]` when `instance` is non-empty but incomprehensive, redundant, non-positive, and mixed, with `blindness = 2_IK`, with `sorted = .false._LK, positive = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
#endif
            pattern = array
            instance = [10_IK, 3_IK]
            loc_ref = [integer(IK)::]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.`.", int(__LINE__, IK))

            instance = [3_IK, 10_IK]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK, `sorted = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK, `sorted = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK, `positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK, `positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK, `sorted = .true._LK, positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK, `sorted = .true._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 2_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` the same as the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 2_IK, `sorted = .false._LK, positive = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            pattern = array//array
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            pattern = [array, array]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            pattern = [array, array]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            pattern = [array, array]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            pattern = [array, array]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            pattern = [array, array]
#endif
            instance = [1_IK, 3_IK]
            loc_ref = [integer(IK)::]

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .true._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `sorted = .false._LK, positive = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 20_IK`, `sorted = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 20_IK`, `sorted = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .true._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 20_IK`, `positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, positive = .false._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 20_IK`, `positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .false._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 20_IK`, `sorted = .true._LK, positive = .false.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .true._LK, positive = .true._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 20_IK`, `sorted = .true._LK, positive = .true.`.", int(__LINE__, IK))

            call report(__LINE__, iseq, instance = instance, sorted = .false._LK, positive = .false._LK, blindness = 20_IK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `pattern` larger than the `array` must yield `loc = [integer(IK)::]` when `instance` is ordered, positive, non-empty but non-sensible, with `blindness = 20_IK`, `sorted = .false._LK, positive = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report   ( line & ! LCOV_EXCL_LINE
                            , iseq & ! LCOV_EXCL_LINE
                            , instance & ! LCOV_EXCL_LINE
                            , blindness & ! LCOV_EXCL_LINE
                            , positive & ! LCOV_EXCL_LINE
                            , sorted & ! LCOV_EXCL_LINE
                            )
            integer, intent(in) :: line
            logical(LK) , external  , optional              :: iseq
            integer(IK) , intent(in), optional, contiguous  :: instance(:)
            integer(IK) , intent(in), optional              :: blindness
            logical(LK) , intent(in), optional              :: positive
            logical(LK) , intent(in), optional              :: sorted
            logical(LK) :: sorted_def, positive_def
            integer(IK) :: blindness_def, nloc

            sorted_def = .false._LK; if (present(sorted)) sorted_def = sorted
            positive_def = .false._LK; if (present(positive)) positive_def = positive
            blindness_def = 1_IK; if (present(blindness)) blindness_def = blindness

            if (setLocEnabled) then
                if (present(iseq) .and. present(instance)) then
                    call setLoc(loc, nloc, array, pattern, iseq, instance, sorted_def, positive_def, blindness_def)
                elseif (present(instance)) then
                    call setLoc(loc, nloc, array, pattern, instance, sorted_def, positive_def, blindness_def)
                elseif (present(iseq)) then
                    call setLoc(loc, nloc, array, pattern, iseq, blindness_def)
                else
                    call setLoc(loc, nloc, array, pattern, blindness_def)
                end if
            else
                if (present(iseq) .and. present(instance) .and. present(sorted) .and. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, instance = instance, sorted = sorted, positive = positive, blindness = blindness)
                elseif (present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, instance = instance, sorted = sorted, blindness = blindness)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, instance = instance, positive = positive, blindness = blindness)
                elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, instance = instance, blindness = blindness)
                elseif (present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, iseq = iseq, blindness = blindness)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. present(positive)) then
                    loc = getLoc(array, pattern, instance = instance, sorted = sorted, positive = positive, blindness = blindness)
                elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, instance = instance, sorted = sorted, blindness = blindness)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(positive)) then
                    loc = getLoc(array, pattern, instance = instance, positive = positive, blindness = blindness)
                elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, instance = instance, blindness = blindness)
                elseif (.not. present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(positive)) then
                    loc = getLoc(array, pattern, blindness = blindness)
                else
                    error stop PROCEDURE_NAME//SK_": Unrecognized interface in testing." ! LCOV_EXCL_LINE
                end if
                nloc = size(loc, 1, IK)
            end if

            ! Report test results if needed.

            !write(*,*) setLocEnabled, present(instance), present(sorted), present(positive)
            !write(*,*) array
            !write(*,*) pattern
            !write(*,*) loc
            !write(*,*) loc_ref

            assertion = assertion .and. nloc == size(loc_ref, 1, IK)
            if (assertion .and. 0_IK < nloc) assertion = assertion .and. all(loc(1 : nloc) == loc_ref)

            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("setLocEnabled")
                call test%disp%show( setLocEnabled )
                call test%disp%show("array")
                call test%disp%show( array )
                call test%disp%show("pattern")
                call test%disp%show( pattern )
                call test%disp%show("nloc")
                call test%disp%show( nloc )
                call test%disp%show("loc_ref")
                call test%disp%show( loc_ref )
                call test%disp%show("loc(1:nloc)")
                call test%disp%show( loc(1:nloc) )
                call test%disp%show("present(iseq)")
                call test%disp%show( present(iseq) )
                call test%disp%show("present(instance)")
                call test%disp%show( present(instance) )
                if (present(instance)) then
                    call test%disp%show("instance")
                    call test%disp%show( instance )
                end if
                call test%disp%show("present(positive)")
                call test%disp%show( present(positive) )
                if (present(positive)) then
                    call test%disp%show("positive")
                    call test%disp%show( positive )
                end if
                call test%disp%show("present(sorted)")
                call test%disp%show( present(sorted) )
                if (present(sorted)) then
                    call test%disp%show("sorted")
                    call test%disp%show( sorted )
                end if
                call test%disp%show("present(blindness)")
                call test%disp%show( present(blindness) )
                if (present(blindness)) then
                    call test%disp%show("blindness")
                    call test%disp%show( blindness )
                end if
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if

            if (nloc /= size(loc_ref, 1, IK)) call test%assert(assertion, PROCEDURE_NAME//SK_": The condition `nloc = size(loc_ref) must hold.", int(line, IK))

        end subroutine
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef IS_EQUAL
#undef ALL