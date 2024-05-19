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
!>  This module contains implementations of the tests of the procedures of [pm_arrayMerge](@ref pm_arrayMerge).
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMerged_D0_SK_ENABLED || setMerged_D0_SK_ENABLED
#define GET_SIZE(array) len(array, kind = IK)
#define GET_INDEX(i) i:i
#else
#define GET_SIZE(array) size(array, kind = IK)
#define GET_INDEX(i) i
#endif

#if     getMerged_D1_LK_ENABLED || setMerged_D1_LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     getMerged_D1_PSSK_ENABLED || setMerged_D1_PSSK_ENABLED
        use pm_container, only: strc => css_pdt, operator(==)
#endif
        integer(IK) :: i

#if     getMerged_D0_SK_ENABLED || setMerged_D0_SK_ENABLED
        character(:,SKG), allocatable   :: SortedArray1, SortedArray2, SortedMergedArray, SortedMergedArray_ref
#elif   getMerged_D1_SK_ENABLED || setMerged_D1_SK_ENABLED
        character(2,SKG), allocatable   :: SortedArray1(:), SortedArray2(:), SortedMergedArray(:), SortedMergedArray_ref(:)
#elif   getMerged_D1_IK_ENABLED || setMerged_D1_IK_ENABLED
        integer(IKG)    , allocatable   :: SortedArray1(:), SortedArray2(:), SortedMergedArray(:), SortedMergedArray_ref(:)
#elif   getMerged_D1_LK_ENABLED || setMerged_D1_LK_ENABLED
        logical(LKG)    , allocatable   :: SortedArray1(:), SortedArray2(:), SortedMergedArray(:), SortedMergedArray_ref(:)
#elif   getMerged_D1_CK_ENABLED || setMerged_D1_CK_ENABLED
        complex(CKG)    , allocatable   :: SortedArray1(:), SortedArray2(:), SortedMergedArray(:), SortedMergedArray_ref(:)
#elif   getMerged_D1_RK_ENABLED || setMerged_D1_RK_ENABLED
        real(RKG)       , allocatable   :: SortedArray1(:), SortedArray2(:), SortedMergedArray(:), SortedMergedArray_ref(:)
#elif   getMerged_D1_PSSK_ENABLED || setMerged_D1_PSSK_ENABLED
        type(strc(SKG)) , allocatable   :: SortedArray1(:), SortedArray2(:), SortedMergedArray(:), SortedMergedArray_ref(:)
#else
#error  "Unrecognized Interface."
#endif

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getMerged_D1_IK_ENABLED || setMerged_D1_IK_ENABLED || \
        getMerged_D1_RK_ENABLED || setMerged_D1_RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     getMerged_D1_IK_ENABLED || setMerged_D1_IK_ENABLED
        SortedArray1 = int([1, 2], kind = IKG)
        SortedArray2 = int([-1], kind = IKG)
        SortedMergedArray_ref = int([-1, 1, 2], IKG)
#elif   getMerged_D1_RK_ENABLED || setMerged_D1_RK_ENABLED
        SortedArray1 = [1._RKG, 2._RKG]
        SortedArray2 = [-1._RKG]
        SortedMergedArray_ref = [-1._RKG, 1._RKG, 2._RKG]
#endif
        call report(int(__LINE__, IK))
        call test%assert(assertion, SK_"The procedure must be able to merge-sort two input contiguous sorted arrays of rank 1 where all elements of array2 are smaller than elements of array1.", int(__LINE__, IK))
        call test%assert(assertion, SK_"The procedure must be able to merge-sort two input contiguous sorted arrays of rank 1 where all elements of array2 are smaller than elements of array1.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     getMerged_D1_IK_ENABLED || setMerged_D1_IK_ENABLED
        SortedArray1 = [1_IKG, 4_IKG]
        SortedArray2 = [2_IKG, 3_IKG]
        SortedMergedArray_ref = [1_IKG, 2_IKG, 3_IKG, 4_IKG]
#elif   getMerged_D1_RK_ENABLED || setMerged_D1_RK_ENABLED
        SortedArray1 = [1._RKG, 4._RKG]
        SortedArray2 = [2._RKG, 3._RKG]
        SortedMergedArray_ref = [1._RKG, 2._RKG, 3._RKG, 4._RKG]
#endif
        call report(int(__LINE__, IK))
        call test%assert(assertion, SK_"The procedure must be able to merge-sort two input contiguous sorted arrays of rank 1 where all elements of array2 are between elements of array1.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     getMerged_D1_IK_ENABLED || setMerged_D1_IK_ENABLED
        SortedArray1 = [1_IKG, 4_IKG, 10_IKG]
        SortedArray2 = [-2_IKG, 3_IKG, 7_IKG, 11_IKG]
        SortedMergedArray_ref = [-2_IKG, 1_IKG, 3_IKG, 4_IKG, 7_IKG, 10_IKG, 11_IKG]
#elif   getMerged_D1_RK_ENABLED || setMerged_D1_RK_ENABLED
        SortedArray1 = [1._RKG, 4._RKG, 10._RKG]
        SortedArray2 = [-2._RKG, 3._RKG, 7._RKG, 11._RKG]
        SortedMergedArray_ref = [-2._RKG, 1._RKG, 3._RKG, 4._RKG, 7._RKG, 10._RKG, 11._RKG]
#endif
        call report(int(__LINE__, IK))
        call test%assert(assertion, SK_"The procedure must be able to merge-sort two input arbitrarily-sized and valued contiguous sorted arrays of rank 1.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

        if (allocated(SortedArray1)) deallocate(SortedArray1) ! LCOV_EXCL_LINE
        allocate(SortedArray1(0))
#if     getMerged_D1_IK_ENABLED || setMerged_D1_IK_ENABLED
        SortedArray2 = [-2_IKG, 3_IKG, 7_IKG, 11_IKG]
        SortedMergedArray_ref = SortedArray2
#elif   getMerged_D1_RK_ENABLED || setMerged_D1_RK_ENABLED
        SortedArray2 = [-2._RKG, 3._RKG, 7._RKG, 11._RKG]
        SortedMergedArray_ref = SortedArray2
#endif
        call report(int(__LINE__, IK))
        call test%assert(assertion, SK_"The procedure must be able to merge-sort an empty array with a non-empty array.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

        if (allocated(SortedArray2)) deallocate(SortedArray2) ! LCOV_EXCL_LINE
        allocate(SortedArray2(0))
#if     getMerged_D1_IK_ENABLED || setMerged_D1_IK_ENABLED
        SortedArray1 = [-2_IKG, 3_IKG, 7_IKG, 11_IKG]
        SortedMergedArray_ref = SortedArray1
#elif   getMerged_D1_RK_ENABLED || setMerged_D1_RK_ENABLED
        SortedArray1 = [-2._RKG, 3._RKG, 7._RKG, 11._RKG]
        SortedMergedArray_ref = SortedArray1
#endif
        call report(int(__LINE__, IK))
        call test%assert(assertion, SK_"The procedure must be able to merge-sort a non-empty array with an empty array.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

        if (allocated(SortedArray2)) deallocate(SortedArray2) ! LCOV_EXCL_LINE
        allocate(SortedArray2(0))
#if     getMerged_D1_IK_ENABLED || setMerged_D1_IK_ENABLED
        SortedArray1 = [-2_IKG, 3_IKG, 7_IKG, 11_IKG]
        SortedMergedArray_ref = SortedArray1
#elif   getMerged_D1_RK_ENABLED || setMerged_D1_RK_ENABLED
        SortedArray1 = [-2._RKG, 3._RKG, 7._RKG, 11._RKG]
        SortedMergedArray_ref = SortedArray1
#endif
        call report(int(__LINE__, IK))
        call test%assert(assertion, SK_"The procedure must be able to merge-sort a non-empty array with an empty array.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call testGeneric()
        call testGeneric(isSorted)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine testGeneric(isSorted)
            procedure(logical(LK))  , optional  :: isSorted
            do i = 1, 200
                call reset()
#if             getMerged_D0_SK_ENABLED || setMerged_D0_SK_ENABLED
                allocate(character(getUnifRand(0, 10),SKG) :: SortedArray1)
                allocate(character(getUnifRand(0, 10),SKG) :: SortedArray2)
                call setUnifRand(SortedArray1)
                call setUnifRand(SortedArray2)
                SortedMergedArray_ref = SortedArray1//SortedArray2
#elif           getMerged_D1_SK_ENABLED || setMerged_D1_SK_ENABLED
                allocate(SortedArray1(getUnifRand(0, 10)))
                allocate(SortedArray2(getUnifRand(0, 10)))
                allocate(SortedMergedArray_ref(size(SortedArray1) + size(SortedArray2)))
                call setUnifRand(SortedArray1)
                call setUnifRand(SortedArray2)
                SortedMergedArray_ref(1:size(SortedArray1)) = SortedArray1
                SortedMergedArray_ref(size(SortedArray1) + 1 : size(SortedArray1) + size(SortedArray2)) = SortedArray2
#elif           getMerged_D1_PSSK_ENABLED || setMerged_D1_PSSK_ENABLED
                allocate(SortedArray1(getUnifRand(0, 10)))
                allocate(SortedArray2(getUnifRand(0, 10)))
                allocate(SortedMergedArray_ref(size(SortedArray1) + size(SortedArray2)))
                block
                    integer :: j
                    do j = 1, size(SortedArray1, kind = IK)
                        call setUnifRand(SortedArray1(j)%val)
                    end do
                    do j = 1, size(SortedArray1, kind = IK)
                        call setUnifRand(SortedArray2(j)%val)
                    end do
                    SortedMergedArray_ref(1:size(SortedArray1)) = SortedArray1
                    SortedMergedArray_ref(size(SortedArray1) + 1 : size(SortedArray1) + size(SortedArray2)) = SortedArray2
                end block
#else
                allocate(SortedArray1(getUnifRand(0, 10)))
                allocate(SortedArray2(getUnifRand(0, 10)))
                call setUnifRand(SortedArray1)
                call setUnifRand(SortedArray2)
                SortedMergedArray_ref = [SortedArray1, SortedArray2]
#endif
                if (present(isSorted)) then
                    call setSorted(SortedMergedArray_ref, isSorted)
                    call setSorted(SortedArray1, isSorted)
                    call setSorted(SortedArray2, isSorted)
                else
                    call setSorted(SortedMergedArray_ref)
                    call setSorted(SortedArray1)
                    call setSorted(SortedArray2)
                end if
                call report(int(__LINE__, IK), isSorted)
                call test%assert(assertion, SK_"The procedure must be able to merge-sort two arrays to return a merged array with present(isSorted) = "//getStr(present(isSorted))//SK_".", int(__LINE__, IK))
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function isSorted(lhs, rhs) result(sorted)
#if         getMerged_D0_SK_ENABLED || setMerged_D0_SK_ENABLED
            character(1,SKG)        , intent(in) :: lhs, rhs
#elif       getMerged_D1_SK_ENABLED || setMerged_D1_SK_ENABLED
            character(*,SKG)        , intent(in) :: lhs, rhs
#elif       getMerged_D1_IK_ENABLED || setMerged_D1_IK_ENABLED
            integer(IKG)            , intent(in) :: lhs, rhs
#elif       getMerged_D1_LK_ENABLED || setMerged_D1_LK_ENABLED
            use pm_logicalCompare, only: operator(>)
            logical(LKG)            , intent(in) :: lhs, rhs
#elif       getMerged_D1_CK_ENABLED || setMerged_D1_CK_ENABLED
            use pm_complexCompareLex, only: operator(>)
            complex(CKG)            , intent(in) :: lhs, rhs
#elif       getMerged_D1_RK_ENABLED || setMerged_D1_RK_ENABLED
            real(RKG)               , intent(in) :: lhs, rhs
#elif       getMerged_D1_PSSK_ENABLED || setMerged_D1_PSSK_ENABLED
            use pm_container, only: css_pdt, operator(>)
            type(css_pdt(SKG)) , intent(in) :: lhs, rhs
#else
#error      "Unrecognized interface."
#endif
            logical(LK) :: sorted
            sorted = lhs > rhs
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(SortedArray1)) deallocate(SortedArray1)
            if (allocated(SortedArray2)) deallocate(SortedArray2)
            if (allocated(SortedMergedArray)) deallocate(SortedMergedArray)
            if (allocated(SortedMergedArray_ref)) deallocate(SortedMergedArray_ref)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, isSorted)
            integer(IK), intent(in) :: line
            procedure(logical(LK))  , optional  :: isSorted
            integer(IK) :: i
#if         getMerged_ENABLED
            if (present(isSorted)) then
                SortedMergedArray = getMerged(SortedArray1, SortedArray2, isSorted)
            else
                SortedMergedArray = getMerged(SortedArray1, SortedArray2)
            end if
#elif       setMerged_ENABLED
            if (allocated(SortedMergedArray)) deallocate(SortedMergedArray)
#if         getMerged_D0_SK_ENABLED || setMerged_D0_SK_ENABLED
            allocate(character(len(SortedMergedArray_ref),SKG) :: SortedMergedArray)
#else
            allocate(SortedMergedArray, mold = SortedMergedArray_ref)
#endif
            if (present(isSorted)) then
                call setMerged(SortedMergedArray, SortedArray1, SortedArray2, isSorted)
            else
                call setMerged(SortedMergedArray, SortedArray1, SortedArray2)
            end if
#else
#error      "Unrecognized interface."
#endif
            assertion = assertion .and. GET_SIZE(SortedMergedArray) == GET_SIZE(SortedMergedArray_ref)
            call test%assert(assertion, SK_"The sizes of `SortedMergedArray` and `SortedMergedArray_ref` must equal with present(isSorted) = "//getStr(present(isSorted))//SK_": "//getStr([GET_SIZE(SortedMergedArray), GET_SIZE(SortedMergedArray_ref)]), line)

            do i = 1, GET_SIZE(SortedMergedArray)
#if             getMerged_D1_LK_ENABLED || setMerged_D1_LK_ENABLED
                assertion = assertion .and. (logical(SortedMergedArray(GET_INDEX(i))) IS_EQUAL logical(SortedMergedArray_ref(GET_INDEX(i))))
#else
                assertion = assertion .and. (SortedMergedArray(GET_INDEX(i)) IS_EQUAL SortedMergedArray_ref(GET_INDEX(i)))
#endif
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0,:,', '))")
                    write(test%disp%unit,"(*(g0,:,', '))") "i", i
                    write(test%disp%unit,"(*(g0,:,', '))") "SortedMergedArray(i)", SortedMergedArray(GET_INDEX(i))
                    write(test%disp%unit,"(*(g0,:,', '))") "SortedMergedArray_ref(i)", SortedMergedArray_ref(GET_INDEX(i))
                    write(test%disp%unit,"(*(g0,:,', '))") "SortedMergedArray(GET_INDEX(i)) IS_EQUAL SortedMergedArray_ref(GET_INDEX(i))", SortedMergedArray(GET_INDEX(i)) IS_EQUAL SortedMergedArray_ref(GET_INDEX(i))
#if                 getMerged_D1_LK_ENABLED || setMerged_D1_LK_ENABLED
                    write(test%disp%unit,"(*(g0,:,', '))") SortedMergedArray
                    write(test%disp%unit,"(*(g0,:,', '))") SortedMergedArray_ref
                    write(test%disp%unit,"(*(g0,:,', '))") SortedMergedArray(GET_INDEX(i)) .eqv. SortedMergedArray_ref(GET_INDEX(i))
#endif
                    write(test%disp%unit,"(*(g0,:,', '))")
                    return
                    ! LCOV_EXCL_STOP
                end if
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  GET_INDEX
#undef  IS_EQUAL
#undef  GET_SIZE