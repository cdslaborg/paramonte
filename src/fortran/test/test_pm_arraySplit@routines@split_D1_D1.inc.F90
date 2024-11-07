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
!>  This module contains implementations of the tests of the procedures under the generic interfaces of [pm_arraySplit](@ref pm_arraySplit).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED && D1_D1_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     SK_ENABLED && D0_D0_ENABLED
#define GET_SIZE(array) len(array, kind = IK)
#define ALL
#else
#define GET_SIZE(array) size(array, kind = IK)
#endif
        character(*, SK), parameter :: PROCEDURE_NAME = "@split()"
#if     SK_ENABLED && D0_D0_ENABLED
#if     CVXK_ENABLED
        type(css_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(css_pdt(SKG)), allocatable :: field(:), field_ref(:)
#endif
        character(:,SKG), allocatable :: array, sep
#elif   SK_ENABLED && D1_D1_ENABLED
#if     CVXK_ENABLED
        type(cvs_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvs_pdt(SKG)), allocatable :: field(:), field_ref(:)
#endif
        character(2,SKG), allocatable :: array(:), sep(:)
#elif   IK_ENABLED && D1_D1_ENABLED
#if     CVXK_ENABLED
        type(cvi_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvi_pdt(IKG)), allocatable :: field(:), field_ref(:)
#endif
        integer(IKG), allocatable :: array(:), sep(:)
#elif   LK_ENABLED && D1_D1_ENABLED
#if     CVXK_ENABLED
        type(cvl_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvl_pdt(LKG)), allocatable :: field(:), field_ref(:)
#endif
        logical(LKG), allocatable :: array(:), sep(:)
#elif   CK_ENABLED && D1_D1_ENABLED
#if     CVXK_ENABLED
        type(cvc_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvc_pdt(CKG)), allocatable :: field(:), field_ref(:)
#endif
        complex(CKG), allocatable :: array(:), sep(:)
#elif   RK_ENABLED && D1_D1_ENABLED
#if     CVXK_ENABLED
        type(cvr_type), allocatable :: field(:), field_ref(:)
#elif   PVXK_ENABLED
        type(cvr_pdt(RKG)), allocatable :: field(:), field_ref(:)
#endif
        real(RKG), allocatable :: array(:), sep(:)
#else
#error  "Unrecognized interface."
#endif
        integer(IK), allocatable :: sindex(:,:), sindex_ref(:,:)
        integer(IK), allocatable :: instance(:)
        integer(IK) :: i

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        do i = 1, 2 ! repeat for index splitting.
            call runTestsWith()
            call runTestsWith(iseq = iseq)
            call runTestsWith(keep = .true._LK)
            call runTestsWith(keep = .false._LK)
            call runTestsWith(iseq = iseq, keep = .true._LK)
            call runTestsWith(iseq = iseq, keep = .false._LK)
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function iseq(Segment, sep & !LCOV_EXCL_LINE
#if     !(SK_ENABLED && D0_D0_ENABLED)
        , lenDelim & !LCOV_EXCL_LINE
#endif
        ) result(equivalent)
            logical(LK) :: equivalent
#if         !(SK_ENABLED && D0_D0_ENABLED)
            integer(IK), intent(in) :: lenDelim
#endif
#if         SK_ENABLED && D0_D0_ENABLED
            character(*,SKG), intent(in) :: sep, Segment
#elif       SK_ENABLED && D1_D1_ENABLED
            character(*,SKG), intent(in) :: Segment(lenDelim), sep(lenDelim)
#elif       IK_ENABLED && D1_D1_ENABLED
            integer(IKG)    , intent(in) :: Segment(lenDelim), sep(lenDelim)
#elif       LK_ENABLED && D1_D1_ENABLED
            logical(LKG)    , intent(in) :: Segment(lenDelim), sep(lenDelim)
#elif       CK_ENABLED && D1_D1_ENABLED
            complex(CKG)    , intent(in) :: Segment(lenDelim), sep(lenDelim)
#elif       RK_ENABLED && D1_D1_ENABLED
            real(RKG)       , intent(in) :: Segment(lenDelim), sep(lenDelim)
#endif
            equivalent = ALL(sep IS_EQUAL Segment)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(array)) deallocate(array)
            if (allocated(instance)) deallocate(instance)
            if (allocated(sep)) deallocate(sep)
            if (allocated(sindex)) deallocate(sindex)
            if (allocated(sindex_ref)) deallocate(sindex_ref)
#if         CVXK_ENABLED || PVXK_ENABLED
            if (allocated(field_ref)) deallocate(field_ref) ! \bug gfortran 11 will not even compile this line.
            if (allocated(field)) deallocate(field) ! \bug gfortran 11 will not even compile this line.
#endif
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(iseq, keep)
            logical(LK) , intent(in), optional  :: keep
            logical(LK) , external  , optional  :: iseq

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = ""
            sep = " "
#elif       SK_ENABLED && D1_D1_ENABLED
            sep = [" "]
            allocate(array(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            sep = [1_IKG]
            allocate(array(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            sep = [.false._LK]
            allocate(array(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            sep = [1._CKG]
            allocate(array(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            sep = [1._RKG]
            allocate(array(0))
#endif
#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif

            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1
            sindex_ref(2,1) = 0

            allocate(instance(0))

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty array, has empty resulting array with `instance` with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = ""
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            allocate(character(2) :: sep(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            allocate(sep(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            allocate(sep(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [1._CKG, 1._CKG]
            allocate(sep(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            allocate(sep(0))
#endif
#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1
            sindex_ref(2,1) = 2

            allocate(instance(0))

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is empty.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is empty with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is empty with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is empty with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is empty with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is empty with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is empty with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is empty with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is empty with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = ""
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            allocate(character(2) :: sep(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            allocate(sep(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            allocate(sep(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [1._CKG, 1._CKG]
            allocate(sep(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            allocate(sep(0))
#endif
#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1
            sindex_ref(2,1) = 2

            instance = [1_IK]

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is non-empty.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is non-empty with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is non-empty with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is non-empty with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is non-empty with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is non-empty with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is non-empty with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is non-empty with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with empty `sep` must yield an `field` of size 1 whose sole component is `array`-valued when `instance` is non-empty with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            endif
#endif
#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "A"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = ""
                field_ref(2)%val = sep
                field_ref(3)%val = ""
                field_ref(4)%val = sep
                field_ref(5)%val = ""
            else
                field_ref(1)%val = ""
                field_ref(2)%val = ""
                field_ref(3)%val = ""
            endif
#endif
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["AA"]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(character(2,SKG) :: field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(character(2,SKG) :: field_ref(5)%val(0))
            else
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                allocate(character(2,SKG) :: field_ref(2)%val(0))
                allocate(character(2,SKG) :: field_ref(3)%val(0))
            endif
#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [1_IKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.false._LK]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG)]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [1._RKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
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

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            endif
#endif
#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "A"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = ""
                field_ref(2)%val = sep
                field_ref(3)%val = ""
                field_ref(4)%val = sep
                field_ref(5)%val = ""
            else
                field_ref(1)%val = ""
                field_ref(2)%val = ""
                field_ref(3)%val = ""
            endif
#endif
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["AA"]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(character(2,SKG) :: field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(character(2,SKG) :: field_ref(5)%val(0))
            else
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                allocate(character(2,SKG) :: field_ref(2)%val(0))
                allocate(character(2,SKG) :: field_ref(3)%val(0))
            endif
#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [1_IKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.false._LK]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG)]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [1._RKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
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

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and redundant with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            endif
#endif
#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "A"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = ""
                field_ref(2)%val = sep
                field_ref(3)%val = ""
                field_ref(4)%val = sep
                field_ref(5)%val = ""
            else
                field_ref(1)%val = ""
                field_ref(2)%val = ""
                field_ref(3)%val = ""
            endif
#endif
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["AA"]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(character(2,SKG) :: field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(character(2,SKG) :: field_ref(5)%val(0))
            else
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                allocate(character(2,SKG) :: field_ref(2)%val(0))
                allocate(character(2,SKG) :: field_ref(3)%val(0))
            endif
#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [1_IKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.false._LK]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG)]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [1._RKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
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

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .false._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty and mixed with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            endif
#endif
#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "A"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = ""
                field_ref(2)%val = sep
                field_ref(3)%val = ""
                field_ref(4)%val = sep
                field_ref(5)%val = ""
            else
                field_ref(1)%val = ""
                field_ref(2)%val = ""
                field_ref(3)%val = ""
            endif
#endif
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["AA"]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(character(2,SKG) :: field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(character(2,SKG) :: field_ref(5)%val(0))
            else
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                allocate(character(2,SKG) :: field_ref(2)%val(0))
                allocate(character(2,SKG) :: field_ref(3)%val(0))
            endif
#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [1_IKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.false._LK]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG)]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [1._RKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
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

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
            else
                allocate(field_ref(3))
            endif
#endif
#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "A"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = ""
                field_ref(2)%val = sep
                field_ref(3)%val = ""
                field_ref(4)%val = sep
                field_ref(5)%val = ""
            else
                field_ref(1)%val = ""
                field_ref(2)%val = ""
                field_ref(3)%val = ""
            endif
#endif
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["AA"]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(character(2,SKG) :: field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(character(2,SKG) :: field_ref(5)%val(0))
            else
                allocate(character(2,SKG) :: field_ref(1)%val(0))
                allocate(character(2,SKG) :: field_ref(2)%val(0))
                allocate(character(2,SKG) :: field_ref(3)%val(0))
            endif
#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [1_IKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.false._LK]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG)]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [1._RKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(1)%val(0))
                field_ref(2)%val = sep
                allocate(field_ref(3)%val(0))
                field_ref(4)%val = sep
                allocate(field_ref(5)%val(0))
            else
                allocate(field_ref(1)%val(0))
                allocate(field_ref(2)%val(0))
                allocate(field_ref(3)%val(0))
            endif
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

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, unsorted, and mixed.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size `size(array)/size(sep) + 1` whose components are all empty when `instance` is non-empty, redundant, non-unique, unsorted, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(3))
            else
                allocate(field_ref(2))
            endif
#endif
#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "A"
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = ""
                field_ref(2)%val = sep
                field_ref(3)%val = "A"
            else
                field_ref(1)%val = ""
                field_ref(2)%val = "A"
            endif
#endif
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["AA"]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [character(2,SKG) ::]
                field_ref(2)%val = sep
                field_ref(3)%val = ["AA"]
            else
                field_ref(1)%val = [character(2,SKG) ::]
                field_ref(2)%val = ["AA"]
            endif
#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [1_IKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = sep
                field_ref(3)%val = [1_IKG]
            else
                field_ref(1)%val = [integer(IKG) ::]
                field_ref(2)%val = [1_IKG]
            endif
#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.false._LK]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = sep
                field_ref(3)%val = [.false._LKG]
            else
                field_ref(1)%val = [logical(LKG) ::]
                field_ref(2)%val = [.false._LKG]
            endif
#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG)]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = sep
                field_ref(3)%val = [(1._CKG,-1._CKG)]
            else
                field_ref(1)%val = [complex(CKG) ::]
                field_ref(2)%val = [(1._CKG,-1._CKG)]
            endif
#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [1._RKG]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = sep
                field_ref(3)%val = [1._RKG]
            else
                field_ref(1)%val = [real(RKG) ::]
                field_ref(2)%val = [1._RKG]
            endif
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

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "A"
!#if         CVXK_ENABLED || PVXK_ENABLED
!            field_ref(1)%val = ""
!            field_ref(2)%val = "A"
!#endif
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["AA"]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(character(2,SKG) :: field_ref(1)%val(0))
!            field_ref(2)%val = ["AA"]
!#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [1_IKG]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            field_ref(2)%val = [1_IKG]
!#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.false._LK]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            field_ref(2)%val = [.false._LK]
!#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG)]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            field_ref(2)%val = [(1._CKG,-1._CKG)]
!#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [1._RKG]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            field_ref(2)%val = [1._RKG]
!#endif
#endif
            instance = [-2_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(3))
                field_ref(1)%val = array(1:0)
                field_ref(2)%val = array(1:1)
                field_ref(3)%val = array(2:2)
            else
                allocate(field_ref(2))
                field_ref(1)%val = array(1:0)
                field_ref(2)%val = array(2:2)
            endif
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
            endif

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose first component is empty but the last corresponds to the last element of `array` when `instance` is non-empty but incomprehensive and unique, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "A"
!#if         CVXK_ENABLED || PVXK_ENABLED
!            field_ref(1)%val = "A"
!            field_ref(2)%val = ""
!#endif
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["AA"]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            field_ref(1)%val = ["AA"]
!            allocate(character(2,SKG) :: field_ref(2)%val(0))
!#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [1_IKG]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            field_ref(1)%val = [1_IKG]
!            allocate(field_ref(2)%val(0))
!#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.false._LK]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            field_ref(1)%val = [.false._LK]
!            allocate(field_ref(2)%val(0))
!#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG)]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            field_ref(1)%val = [(1._CKG,-1._CKG)]
!            allocate(field_ref(2)%val(0))
!#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [1._RKG]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            field_ref(1)%val = [1._RKG]
!            allocate(field_ref(2)%val(0))
!#endif
#endif
            instance = [-1_IK, 3_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(3))
                field_ref(1)%val = array(1:1)
                field_ref(2)%val = array(2:2)
                field_ref(3)%val = array(3:2)
            else
                allocate(field_ref(2))
                field_ref(1)%val = array(1:1)
                field_ref(2)%val = array(3:2)
            endif
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
            endif

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` of non-empty `sep` must yield an `field` of size two whose last component is empty but the first corresponds to the first element of `array` when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
!#if         CVXK_ENABLED || PVXK_ENABLED
!            field_ref(1)%val = ""
!            field_ref(2)%val = ""
!#endif
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(character(2,SKG) :: field_ref(1)%val(0))
!            allocate(character(2,SKG) :: field_ref(2)%val(0))
!#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            allocate(field_ref(2)%val(0))
!#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            allocate(field_ref(2)%val(0))
!#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            allocate(field_ref(2)%val(0))
!#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            allocate(field_ref(2)%val(0))
!#endif
#endif
            sep = array
            instance = [-1_IK, 3_IK]

#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(3))
                field_ref(1)%val = array(1:0)
                field_ref(2)%val = array(1:GET_SIZE(array)) ! fpp
                field_ref(3)%val = array(GET_SIZE(array)+1:GET_SIZE(array)) ! fpp
            else
                allocate(field_ref(2))
                field_ref(1)%val = array(1:0)
                field_ref(2)%val = array(GET_SIZE(array):GET_SIZE(array)-1) ! fpp
            endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = GET_SIZE(array) ! fpp
                sindex_ref(1,3) = GET_SIZE(array); sindex_ref(2,3) = GET_SIZE(array) - 1 ! fpp
            else
                allocate(sindex_ref(2,2))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = GET_SIZE(array); sindex_ref(2,2) = GET_SIZE(array) - 1 ! fpp
            endif

            call reportDelim(iseq, keep, instance = instance)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is non-empty but incomprehensive, redundant, non-unique, and mixed, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
#endif
            sep = array
            instance = [0_IK, 3_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1; sindex_ref(2,1) = GET_SIZE(array) ! fpp

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` the same as the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = array//array
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = [array, array]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [array, array]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [array, array]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [array, array]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [array, array]
#endif
            instance = [0_IK, 3_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1; sindex_ref(2,1) = GET_SIZE(array) ! fpp

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            instance = [1_IK]

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` larger than the `array` must yield an `field` of size two whose components are all empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AABBCCAA"
            sep = "AA"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA", "BB", "BB", "CC", "CC", "AA", "AA"]
            sep = ["AA", "AA"]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(character(2,SKG) :: field_ref(1)%val(0))
!            allocate(character(2,SKG) :: field_ref(3)%val(0))
!            field_ref(2)%val = ["BB", "CC"]
!#endif
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG, 2_IKG, 2_IKG, 3_IKG, 3_IKG, 1_IKG, 1_IKG]
            sep = [1_IKG, 1_IKG]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            field_ref(2)%val = [2_IKG, 2_IKG, 3_IKG, 3_IKG]
!            allocate(field_ref(3)%val(0))
!#endif
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK, .true._LK, .true._LK, .true._LK, .true._LK, .false._LK, .false._LK]
            sep = [.false._LK, .false._LK]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            field_ref(2)%val = [.true._LK, .true._LK, .true._LK, .true._LK]
!            allocate(field_ref(3)%val(0))
!#endif
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG), (2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG), (1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            field_ref(2)%val = [(2._CKG,-2._CKG), (2._CKG,-2._CKG), (3._CKG,-3._CKG), (3._CKG,-3._CKG)]
!            allocate(field_ref(3)%val(0))
!#endif
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG, 2._RKG, 2._RKG, 3._RKG, 3._RKG, 1._RKG, 1._RKG]
            sep = [1._RKG, 1._RKG]
!#if         CVXK_ENABLED || PVXK_ENABLED
!            allocate(field_ref(1)%val(0))
!            field_ref(2)%val = [2._RKG, 2._RKG, 3._RKG, 3._RKG]
!            allocate(field_ref(3)%val(0))
!#endif
#endif
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
                field_ref(1)%val = array(1:0)
                field_ref(2)%val = array(1:2)
                field_ref(3)%val = array(3:6)
                field_ref(4)%val = array(7:8)
                field_ref(5)%val = array(9:8)
            else
                allocate(field_ref(3))
                field_ref(1)%val = array(1:0)
                field_ref(2)%val = array(3:6)
                field_ref(3)%val = array(8:7)
            endif
#endif
            if (getOption(.false._LK, keep)) then
                allocate(sindex_ref(2,5))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 1; sindex_ref(2,2) = 2
                sindex_ref(1,3) = 3; sindex_ref(2,3) = 6
                sindex_ref(1,4) = 7; sindex_ref(2,4) = 8
                sindex_ref(1,5) = 9; sindex_ref(2,5) = 8
            else
                allocate(sindex_ref(2,3))
                sindex_ref(1,1) = 1; sindex_ref(2,1) = 0
                sindex_ref(1,2) = 3; sindex_ref(2,2) = 6
                sindex_ref(1,3) = 9; sindex_ref(2,3) = 8
            endif
            instance = [1_IK, -1_IK]

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with `sep` occurrence at the beginning and the end of the `array` must yield an `field` of size three whose first and last elements are empty when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "X"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["XX"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [0_IKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.true._LK]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(0._CKG,0._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [0._RKG]
#endif
            instance = [0_IK, 3_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1; sindex_ref(2,1) = 2

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            instance = [1_IK]

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with non-matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty, and sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AXAXA"
            sep = "X"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "XX", "AA", "XX", "AA"]
            sep = ["XX"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 0_IKG, 1_IKG, 0_IKG, 1_IKG]
            sep = [0_IKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK]
            sep = [.true._LK]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG), (0._CKG,0._CKG), (1._CKG,-1._CKG)]
            sep = [(0._CKG,0._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 0._RKG, 1._RKG, 0._RKG, 1._RKG]
            sep = [0._RKG]
#endif
            instance = [1_IK, -1_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            if (getOption(.false._LK, keep)) then
                allocate(field_ref(5))
                field_ref(1)%val = array(1:1)
                field_ref(2)%val = array(2:2)
                field_ref(3)%val = array(3:3)
                field_ref(4)%val = array(4:4)
                field_ref(5)%val = array(5:5)
            else
                allocate(field_ref(3))
                field_ref(1)%val = array(1:1)
                field_ref(2)%val = array(3:3)
                field_ref(3)%val = array(5:5)
            endif
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

            call reportDelim(iseq, keep)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with two `sep` occurrences must yield an `field` of size three whose three components are the non-matching elements of `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = "AA"
            sep = "A"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = ["AA", "AA"]
            sep = ["AA"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [1_IKG, 1_IKG]
            sep = [1_IKG]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [.false._LK, .false._LK]
            sep = [.false._LK]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [(1._CKG,-1._CKG), (1._CKG,-1._CKG)]
            sep = [(1._CKG,-1._CKG)]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [1._RKG, 1._RKG]
            sep = [1._RKG]
#endif
            instance = [0_IK, 3_IK]
#if         CVXK_ENABLED || PVXK_ENABLED
            allocate(field_ref(1))
            field_ref(1)%val = array
#endif
            allocate(sindex_ref(2,1))
            sindex_ref(1,1) = 1; sindex_ref(2,1) = 2

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .false.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .true._LK, unique = .true._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .true._LK, unique = .true.`.", int(__LINE__, IK))

            call reportDelim(iseq, keep, instance = instance, sorted = .false._LK, unique = .false._LK)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A non-empty `array` with a matching `sep` must yield an `field` of size one whose sole component is `array` when `instance` is ordered, unique, non-empty but non-sensible, with `sorted = .false._LK, unique = .false.`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reportDelim ( iseq & ! LCOV_EXCL_LINE
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
#if         CVXK_ENABLED || PVXK_ENABLED
            integer(IK)                                     :: i
#endif

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

            ! Report test results if needed.

#if         CVXK_ENABLED || PVXK_ENABLED
            do i = 1, size(field_ref)
                assertion = assertion .and. ALL(field(i)%val IS_EQUAL field_ref(i)%val) ! fpp
            end do
#endif
            assertion = assertion .and. all(sindex == sindex_ref)
            if (test%traceable .and. .not. assertion) then

                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")

                call disp%show("sep")
                call disp%show( sep )
                call disp%show("array")
                call disp%show( array )
                call disp%show("sindex")
                call disp%show( sindex )
                call disp%show("sindex_ref")
                call disp%show( sindex_ref )
#if             CVXK_ENABLED || PVXK_ENABLED
                call disp%show("field")
                call disp%show( field )
                call disp%show("field_ref")
                call disp%show( field_ref )
#endif
                write(test%disp%unit,"(*(g0,:,', '))") "present(instance)  ", present(instance)
                write(test%disp%unit,"(*(g0,:,', '))") "present(unique)    ", present(unique)
                write(test%disp%unit,"(*(g0,:,', '))") "present(sorted)    ", present(sorted)
                write(test%disp%unit,"(*(g0,:,', '))") "present(keep)      ", present(keep)
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

                if (present(keep)) then
                write(test%disp%unit,"(*(g0,:,', '))") "keep               ", keep
                end if

                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP

            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GET_SIZE
#undef IS_EQUAL
#undef ALL