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
!>  [setSplit](@ref pm_arraySplit::setSplit).
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
        character(*, SK), parameter :: PROCEDURE_NAME = "@setSplit()"
#if     SK_ENABLED && D1_D0_ENABLED
#if     CVXK_ENABLED
        type(cvs_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvs_pdt(SKG)), allocatable :: field(:), field_ref(:)
#endif
        character(2,SKG), allocatable :: array(:)
        character(2,SKG) :: sep
#elif   IK_ENABLED && D1_D0_ENABLED
#if     CVXK_ENABLED
        type(cvi_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvi_pdt(IKG)), allocatable :: field(:), field_ref(:)
#endif
        integer(IKG), allocatable :: array(:)
        integer(IKG) :: sep
#elif   LK_ENABLED && D1_D0_ENABLED
#if     CVXK_ENABLED
        type(cvl_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvl_pdt(LKG)), allocatable :: field(:), field_ref(:)
#endif
        logical(LKG), allocatable :: array(:)
        logical(LKG) :: sep
#elif   CK_ENABLED && D1_D0_ENABLED
#if     CVXK_ENABLED
        type(cvc_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvc_pdt(CKG)), allocatable :: field(:), field_ref(:)
#endif
        complex(CKG), allocatable :: array(:)
        complex(CKG) :: sep
#elif   RK_ENABLED && D1_D0_ENABLED
#if     CVXK_ENABLED
        type(cvr_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvr_pdt(RKG)), allocatable :: field(:), field_ref(:)
#endif
        real(RKG), allocatable :: array(:)
        real(RKG) :: sep
#else
#error  "Unrecognized interface."
#endif
        integer(IK), allocatable :: sindex(:,:), sindex_ref(:,:)
        integer(IK), allocatable :: instance(:)
        integer(IK)                                     :: i

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        do i = 1, 2 ! repeat for sindex splitting.
            call runTestsWith()
            call runTestsWith(iseq = iseq)
            call runTestsWith(keep = .true._LK)
            call runTestsWith(keep = .false._LK)
            call runTestsWith(iseq = iseq, keep = .true._LK)
            call runTestsWith(iseq = iseq, keep = .false._LK)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function iseq(segment, sep) result(equivalent)
            logical(LK) :: equivalent
#if         SK_ENABLED && D1_D0_ENABLED
            character(*,SKG), intent(in) :: sep, segment
#elif       IK_ENABLED && D1_D0_ENABLED
            integer(IKG)    , intent(in) :: sep, segment
#elif       LK_ENABLED && D1_D0_ENABLED
            logical(LKG)    , intent(in) :: sep, segment
#elif       CK_ENABLED && D1_D0_ENABLED
            complex(CKG)    , intent(in) :: sep, segment
#elif       RK_ENABLED && D1_D0_ENABLED
            real(RKG)       , intent(in) :: sep, segment
#endif
            equivalent = sep IS_EQUAL segment
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(array)) deallocate(array)
            if (allocated(sindex)) deallocate(sindex)
            if (allocated(instance)) deallocate(instance)
            if (allocated(sindex_ref)) deallocate(sindex_ref)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (allocated(field_ref)) deallocate(field_ref) ! \bug gfortran 11 will not even compile this line.
            if (allocated(field)) deallocate(field) ! \bug gfortran 11 will not even compile this line.
#endif
        end subroutine reset

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(iseq, keep)
            logical(LK) , intent(in), optional  :: keep
            logical(LK) , external  , optional  :: iseq

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            sep = " "
#elif       IK_ENABLED && D1_D0_ENABLED
            sep = 1_IKG
#elif       LK_ENABLED && D1_D0_ENABLED
            sep = .false._LK
#elif       CK_ENABLED && D1_D0_ENABLED
            sep = 1._CKG
#elif       RK_ENABLED && D1_D0_ENABLED
            sep = 1._RKG
#endif
            allocate(instance(0))
            allocate(array(0))
#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1
            sindex_ref(2,1) = 0

            call report(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "AA"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [character(2) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [character(2) ::]
            else
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [character(2) ::]
                field_ref(3)%val = [character(2) ::]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 1_IKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [integer(IKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [integer(IKG) ::]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [integer(IKG) ::]
                field_ref(3)%val = [integer(IKG) ::]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .false._LK
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [logical(LKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [logical(LKG) ::]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [logical(LKG) ::]
                field_ref(3)%val = [logical(LKG) ::]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [complex(CKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [complex(CKG) ::]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [complex(CKG) ::]
                field_ref(3)%val = [complex(CKG) ::]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 1._RKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [real(RKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [real(RKG) ::]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [real(RKG) ::]
                field_ref(3)%val = [real(RKG) ::]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,5))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 2; sindex_ref(2,3) = 1
                sindex_ref(1,4) = 2; sindex_ref(2,4) = 2
                sindex_ref(1,5) = 3; sindex_ref(2,5) = 2
            else
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 3; sindex_ref(2,3) = 2
            end if

            instance = [1_IK, 2_IK]

            call report(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "AA"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [character(2) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [character(2) ::]
            else
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [character(2) ::]
                field_ref(3)%val = [character(2) ::]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 1_IKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [integer(IKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [integer(IKG) ::]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [integer(IKG) ::]
                field_ref(3)%val = [integer(IKG) ::]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .false._LK
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [logical(LKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [logical(LKG) ::]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [logical(LKG) ::]
                field_ref(3)%val = [logical(LKG) ::]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [complex(CKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [complex(CKG) ::]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [complex(CKG) ::]
                field_ref(3)%val = [complex(CKG) ::]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 1._RKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [real(RKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [real(RKG) ::]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [real(RKG) ::]
                field_ref(3)%val = [real(RKG) ::]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,5))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 2; sindex_ref(2,3) = 1
                sindex_ref(1,4) = 2; sindex_ref(2,4) = 2
                sindex_ref(1,5) = 3; sindex_ref(2,5) = 2
            else
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 3; sindex_ref(2,3) = 2
            end if

            instance = [0_IK, 1_IK, 2_IK, 3_IK]

            call report(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "AA"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [character(2) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [character(2) ::]
            else
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [character(2) ::]
                field_ref(3)%val = [character(2) ::]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 1_IKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [integer(IKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [integer(IKG) ::]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [integer(IKG) ::]
                field_ref(3)%val = [integer(IKG) ::]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .false._LK
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [logical(LKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [logical(LKG) ::]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [logical(LKG) ::]
                field_ref(3)%val = [logical(LKG) ::]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [complex(CKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [complex(CKG) ::]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [complex(CKG) ::]
                field_ref(3)%val = [complex(CKG) ::]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 1._RKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [real(RKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [real(RKG) ::]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [real(RKG) ::]
                field_ref(3)%val = [real(RKG) ::]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,5))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 2; sindex_ref(2,3) = 1
                sindex_ref(1,4) = 2; sindex_ref(2,4) = 2
                sindex_ref(1,5) = 3; sindex_ref(2,5) = 2
            else
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 3; sindex_ref(2,3) = 2
            end if

            instance = [1_IK, -1_IK]

            call report(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "AA"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [character(2) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [character(2) ::]
            else
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [character(2) ::]
                field_ref(3)%val = [character(2) ::]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 1_IKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [integer(IKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [integer(IKG) ::]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [integer(IKG) ::]
                field_ref(3)%val = [integer(IKG) ::]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .false._LK
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [logical(LKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [logical(LKG) ::]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [logical(LKG) ::]
                field_ref(3)%val = [logical(LKG) ::]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [complex(CKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [complex(CKG) ::]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [complex(CKG) ::]
                field_ref(3)%val = [complex(CKG) ::]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 1._RKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [real(RKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [real(RKG) ::]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [real(RKG) ::]
                field_ref(3)%val = [real(RKG) ::]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,5))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 2; sindex_ref(2,3) = 1
                sindex_ref(1,4) = 2; sindex_ref(2,4) = 2
                sindex_ref(1,5) = 3; sindex_ref(2,5) = 2
            else
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 3; sindex_ref(2,3) = 2
            end if

            instance = [1_IK, -2_IK, -1_IK]

            call report(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "AA"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [character(2) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [character(2) ::]
            else
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [character(2) ::]
                field_ref(3)%val = [character(2) ::]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 1_IKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [integer(IKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [integer(IKG) ::]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [integer(IKG) ::]
                field_ref(3)%val = [integer(IKG) ::]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .false._LK
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [logical(LKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [logical(LKG) ::]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [logical(LKG) ::]
                field_ref(3)%val = [logical(LKG) ::]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [complex(CKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [complex(CKG) ::]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [complex(CKG) ::]
                field_ref(3)%val = [complex(CKG) ::]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 1._RKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [real(RKG) ::]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [real(RKG) ::]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [real(RKG) ::]
                field_ref(3)%val = [real(RKG) ::]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,5))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 2; sindex_ref(2,3) = 1
                sindex_ref(1,4) = 2; sindex_ref(2,4) = 2
                sindex_ref(1,5) = 3; sindex_ref(2,5) = 2
            else
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 3; sindex_ref(2,3) = 2
            end if

            instance = [1_IK, -1_IK, -2_IK, 0_IK]

            call report(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, unsorted, and mixed.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(3))
            else
                allocate(field_ref(2))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "AA"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = ["AA"]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 1_IKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [1_IKG]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .false._LK
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [.false._LK]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [(1._CKG,-1._CKG)]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 1._RKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [1._RKG]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 2; sindex_ref(2,3) = 2
            else
                allocate(sindex_ref(2,2))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 2
            end if

            instance = [1_IK, -2_IK]

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(3))
            else
                allocate(field_ref(2))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "AA"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = ["AA"]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 1_IKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [1_IKG]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .false._LK
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [.false._LK]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [(1._CKG,-1._CKG)]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 1._RKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [sep]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [1._RKG]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 2; sindex_ref(2,3) = 2
            else
                allocate(sindex_ref(2,2))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 2
            end if

            instance = [-2_IK]

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(3))
            else
                allocate(field_ref(2))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "AA"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = ["AA"]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [character(2) ::]
            else
                field_ref(1)%val = ["AA"]
                field_ref(2)%val = [character(2) ::]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 1_IKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [1_IKG]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [integer(IKG) ::]
            else
                field_ref(1)%val = [1_IKG]
                field_ref(2)%val = [integer(IKG) ::]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .false._LK
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [.false._LK]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [logical(LKG) ::]
            else
                field_ref(1)%val = [.false._LK]
                field_ref(2)%val = [logical(LKG) ::]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [(1._CKG,-1._CKG)]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [complex(CKG) ::]
            else
                field_ref(1)%val = [(1._CKG,-1._CKG)]
                field_ref(2)%val = [complex(CKG) ::]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 1._RKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [1._RKG]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [real(RKG) ::]
            else
                field_ref(1)%val = [1._RKG]
                field_ref(2)%val = [real(RKG) ::]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 1
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 2
                sindex_ref(1,3) = 3; sindex_ref(2,3) = 2
            else
                allocate(sindex_ref(2,2))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 1
                sindex_ref(1,2) = 3; sindex_ref(2,2) = 2
            end if

            instance = [-1_IK, 3_IK]

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of scalar `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(3))
            else
                allocate(field_ref(2))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA"]
            sep = array(1)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [character(2) ::]
            else
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [character(2) ::]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG]
            sep = array(1)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [integer(IKG) ::]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [integer(IKG) ::]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK]
            sep = array(1)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [logical(LKG) ::]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [logical(LKG) ::]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG)]
            sep = array(1)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [complex(CKG) ::]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [complex(CKG) ::]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG]
            sep = array(1)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [real(RKG) ::]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [real(RKG) ::]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 1; sindex_ref(2,3) = 0
            else
                allocate(sindex_ref(2,2))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 0
            end if

            instance = [-1_IK, 3_IK]

            call report(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA"]
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG]
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK]
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG]
#endif
            sep = array(1)
            instance = [0_IK, 3_IK]

#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1; sindex_ref(2,1) = 1

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D1_D0_ENABLED
            allocate(character(2) :: array(0))
            sep = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            allocate(array(0))
            sep = 1_IKG
#elif       LK_ENABLED && D1_D0_ENABLED
            allocate(array(0))
            sep = .false._LK
#elif       CK_ENABLED && D1_D0_ENABLED
            allocate(array(0))
            sep = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            allocate(array(0))
            sep = 1._RKG
#endif
            instance = [0_IK, 3_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1; sindex_ref(2,1) = 0


            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            instance = [1_IK]

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "BB", "BB", "CC", "CC", "AA"]
            sep = "AA"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = ["BB", "BB", "CC", "CC"]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [character(2) ::]
            else
                field_ref(1)%val = [character(2) ::]
                field_ref(2)%val = ["BB", "BB", "CC", "CC"]
                field_ref(3)%val = [character(2) ::]
            end if
#endif
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 2_IKG, 2_IKG, 3_IKG, 3_IKG, 1_IKG]
            sep = 1_IKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [2_IKG, 2_IKG, 3_IKG, 3_IKG]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [integer(IKG) ::]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [2_IKG, 2_IKG, 3_IKG, 3_IKG]
                field_ref(3)%val = [integer(IKG) ::]
            end if
#endif
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .true._LK, .true._LK, .true._LK, .true._LK, .false._LK]
            sep = .false._LK
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [.true._LK, .true._LK, .true._LK, .true._LK]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [logical(LKG) ::]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [.true._LK, .true._LK, .true._LK, .true._LK]
                field_ref(3)%val = [logical(LKG) ::]
            end if
#endif
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [(2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG)]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [complex(CKG) ::]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [(2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG)]
                field_ref(3)%val = [complex(CKG) ::]
            end if
#endif
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 2._RKG, 2._RKG, 3._RKG, 3._RKG, 1._RKG]
            sep = 1._RKG
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [sep]
                field_ref(3)%val = [2._RKG, 2._RKG, 3._RKG, 3._RKG]
                field_ref(4)%val = [sep]
                field_ref(5)%val = [real(RKG) ::]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [2._RKG, 2._RKG, 3._RKG, 3._RKG]
                field_ref(3)%val = [real(RKG) ::]
            end if
#endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,5))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 1
                sindex_ref(1,3) = 2; sindex_ref(2,3) = 5
                sindex_ref(1,4) = 6; sindex_ref(2,4) = 6
                sindex_ref(1,5) = 7; sindex_ref(2,5) = 6
            else
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 5
                sindex_ref(1,3) = 7; sindex_ref(2,3) = 6
            end if

            instance = [1_IK, -1_IK]

            call report(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "XX"
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .true._LK
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 0_IKG
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (0._CKG,0._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 0._RKG
#endif
            instance = [0_IK, 3_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1; sindex_ref(2,1) = 2


            call report(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            instance = [1_IK]

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            end if
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "XX", "AA", "XX", "AA"]
            sep = "XX"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 0_IKG, 1_IKG, 0_IKG, 1_IKG]
            sep = 0_IKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK]
            sep = .true._LK
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG)]
            sep = (0._CKG,0._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 0._RKG, 1._RKG, 0._RKG, 1._RKG]
            sep = 0._RKG
#endif
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = array(1:1)
                field_ref(2)%val = [sep]
                field_ref(3)%val = array(3:3)
                field_ref(4)%val = [sep]
                field_ref(5)%val = array(5:5)
            else
                field_ref(1)%val = array(1:1)
                field_ref(2)%val = array(3:3)
                field_ref(3)%val = array(5:5)
            end if
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,5))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 1
                sindex_ref(1,2) = 2; sindex_ref(2,2) = 2
                sindex_ref(1,3) = 3; sindex_ref(2,3) = 3
                sindex_ref(1,4) = 4; sindex_ref(2,4) = 4
                sindex_ref(1,5) = 5; sindex_ref(2,5) = 5
            else
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 1
                sindex_ref(1,2) = 3; sindex_ref(2,2) = 3
                sindex_ref(1,3) = 5; sindex_ref(2,3) = 5
            end if

            instance = [1_IK, -1_IK]

            call report(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
#endif
#if         SK_ENABLED && D1_D0_ENABLED
            array = ["AA", "AA"]
            sep = "AA"
#elif       IK_ENABLED && D1_D0_ENABLED
            array = [1_IKG, 1_IKG]
            sep = 1_IKG
#elif       LK_ENABLED && D1_D0_ENABLED
            array = [.false._LK, .false._LK]
            sep = .false._LK
#elif       CK_ENABLED && D1_D0_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = (1._CKG,-1._CKG)
#elif       RK_ENABLED && D1_D0_ENABLED
            array = [1._RKG, 1._RKG]
            sep = 1._RKG
#endif
            instance = [0_IK, 3_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1; sindex_ref(2,1) = 2

            call report(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call report(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report   ( iseq & ! LCOV_EXCL_LINE
                            , keep & ! LCOV_EXCL_LINE
                            , instance & ! LCOV_EXCL_LINE
                            , sorted & ! LCOV_EXCL_LINE
                            , unique & ! LCOV_EXCL_LINE
                            )
            use pm_io, only: display_type
            logical(LK) , external  , optional              :: iseq
            logical(LK) , intent(in), optional              :: keep
            integer(IK) , intent(in), optional, contiguous  :: instance(:)
            logical(LK) , intent(in), optional              :: sorted
            logical(LK) , intent(in), optional              :: unique
            type(display_type)                              :: disp

            disp = display_type()

            if (present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                call setSplit(sindex, array, sep, iseq = iseq, instance = instance, sorted = sorted, unique = unique, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, iseq = iseq, instance = instance, sorted = sorted, unique = unique, keep = keep)
#endif
            elseif (present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                call setSplit(sindex, array, sep, iseq = iseq, instance = instance, sorted = sorted, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, iseq = iseq, instance = instance, sorted = sorted, keep = keep)
#endif
            elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                call setSplit(sindex, array, sep, iseq = iseq, instance = instance, unique = unique, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, iseq = iseq, instance = instance, unique = unique, keep = keep)
#endif
            elseif (present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                call setSplit(sindex, array, sep, iseq = iseq, instance = instance, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, iseq = iseq, instance = instance, keep = keep)
#endif
            elseif (present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                call setSplit(sindex, array, sep, iseq = iseq, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, iseq = iseq, keep = keep)
#endif
            elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. present(unique)) then
                call setSplit(sindex, array, sep, instance = instance, sorted = sorted, unique = unique, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, instance = instance, sorted = sorted, unique = unique, keep = keep)
#endif
            elseif (.not. present(iseq) .and. present(instance) .and. present(sorted) .and. .not. present(unique)) then
                call setSplit(sindex, array, sep, instance = instance, sorted = sorted, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, instance = instance, sorted = sorted, keep = keep)
#endif
            elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. present(unique)) then
                call setSplit(sindex, array, sep, instance = instance, unique = unique, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, instance = instance, unique = unique, keep = keep)
#endif
            elseif (.not. present(iseq) .and. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                call setSplit(sindex, array, sep, instance = instance, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, instance = instance, keep = keep)
#endif
            elseif (.not. present(iseq) .and. .not. present(instance) .and. .not. present(sorted) .and. .not. present(unique)) then
                call setSplit(sindex, array, sep, keep = keep)
#if             CVXK_ENABLED || PVXK_ENABLED
                call setSplit(field, array, sep, keep = keep)
#endif
            else
                error stop PROCEDURE_NAME//SK_": Unrecognized interface in testing." ! LCOV_EXCL_LINE
            end if

            !block
            !use pm_io, only: display_type
            !type(display_type) :: disp
            !call disp%show("present(keep)")
            !call disp%show( present(keep) )
            !call disp%show("getOption(.false._LK, keep)")
            !call disp%show( getOption(.false._LK, keep) )
            !call disp%show("field")
            !call disp%show( field )
            !call disp%show("field_ref")
            !call disp%show( field_ref )
            !end block

            ! Report test results if needed.

#if         CVXK_ENABLED || PVXK_ENABLED
            block
                integer(IK) :: isplit
                do isplit = 1_IK, size(field_ref, kind = IK)
                    assertion = assertion .and. ALL(field(isplit)%val IS_EQUAL field_ref(isplit)%val) ! fpp
                end do
            end block
#endif
            assertion = assertion .and. all(sindex == sindex_ref)
            if (test%traceable .and. .not. assertion) then

                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")

                call disp%show("sep")
                call disp%show( sep )
                call disp%show("array")
                call disp%show( array )
#if             CVXK_ENABLED || PVXK_ENABLED
                call disp%show("field")
                call disp%show( field )
                call disp%show("field_ref")
                call disp%show( field_ref )
#endif
                call disp%show("sindex")
                call disp%show( sindex )
                call disp%show("sindex_ref")
                call disp%show( sindex_ref )
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

                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP

            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IS_EQUAL
#undef ALL