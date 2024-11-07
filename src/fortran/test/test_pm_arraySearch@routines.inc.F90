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
!>  [getBin](@ref pm_arraySearch::getBin).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE len
#else
#define GET_SIZE size
#endif

#if     LK_ENABLED && D1_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

        use pm_val2str, only: getStr
        use pm_arraySort, only: setSorted

        character(*, SK), parameter :: PROCEDURE_NAME = "@getBin()"

#if     SK_ENABLED && D0_ENABLED
        character(:,SKG)              , allocatable :: Array
        character(:,SKG)              , allocatable :: value
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: Array
        character(2,SKG)                            :: value
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: Array
        integer(IKG)                                :: value
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: Array
        logical(LKG)                                :: value
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: Array
        complex(CKG)                                :: value
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: Array
        real(RKG)                                   :: value
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: index
        integer(IK) :: index_ref

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK

        call runTestsWith()
        call runTestsWith(isLess = isLess)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function isLess(value, segment) result(less)
            use pm_complexCompareLex, only: operator(<)
#if         SK_ENABLED && D0_ENABLED || SK_ENABLED && D1_ENABLED
            character(*,SKG), intent(in) :: value, segment
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in) :: value, segment
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in) :: value, segment
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in) :: value, segment
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in) :: value, segment
#endif
            logical(LK) :: less
            less = value < segment
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(Array)) deallocate(Array)
        end subroutine reset


        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(isLess)
            logical(LK), external, optional  :: isLess

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_ENABLED
            Array = SKG_"abcdefgh"
            value = SKG_" "
#elif       SK_ENABLED && D1_ENABLED
            Array = ["aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh"]
            value = "  "
#elif       IK_ENABLED && D1_ENABLED
            Array = int([1,2,3,4,5,6,7,8], kind = IKG)
            value = 0_IKG
#elif       CK_ENABLED && D1_ENABLED
            Array = cmplx([1,2,3,4,5,6,7,8], [1,2,3,4,5,6,7,8], kind = CKG)
            value = cmplx(0, kind = CKG)
#elif       RK_ENABLED && D1_ENABLED
            Array = real([1,2,3,4,5,6,7,8], kind = RKG)
            value = 0._RKG
#endif
            index_ref = 0_IK

            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A value that is less than all elements of Array, has `index = 0_IK`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_ENABLED
            Array = "acegikmo"
            value = "a"
#elif       SK_ENABLED && D1_ENABLED
            Array = ["aa", "cc", "ee", "gg", "ii", "kk", "mm", "oo"]
            value = "aa"
#elif       IK_ENABLED && D1_ENABLED
            Array = int([1,2,4,6,8,10,12,14], kind = IKG)
            value = 1_IKG
#elif       CK_ENABLED && D1_ENABLED
            Array = cmplx([1,2,4,6,8,10,12,14], -[1,2,4,6,8,10,12,14], kind = CKG)
            value = cmplx(1, -1, kind = CKG)
#elif       RK_ENABLED && D1_ENABLED
            Array = real([1,2,4,6,8,10,12,14], kind = RKG)
            value = 1._RKG
#endif
            index_ref = 1_IK

            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A value that is equal to the first element of Array, has `index = 1_IK`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_ENABLED
            Array = "acegikmo"
            value = "b"
#elif       SK_ENABLED && D1_ENABLED
            Array = ["aa", "cc", "ee", "gg", "ii", "kk", "mm", "oo"]
            value = "bb"
#elif       IK_ENABLED && D1_ENABLED
            Array = int([0,2,4,6,8,10,12,14], kind = IKG)
            value = 1_IKG
#elif       CK_ENABLED && D1_ENABLED
            Array = cmplx([0,2,4,6,8,10,12,14], -[0,2,4,6,8,10,12,14], kind = CKG)
            value = cmplx(1, 5, kind = CKG)
#elif       RK_ENABLED && D1_ENABLED
            Array = real([0,2,4,6,8,10,12,14], kind = RKG)
            value = 1._RKG
#endif
            index_ref = 1_IK

            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A value that is between the first two elements of Array, has `index = 1_IK`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_ENABLED
            Array = "acegikmo"
            value = "c"
#elif       SK_ENABLED && D1_ENABLED
            Array = ["aa", "cc", "ee", "gg", "ii", "kk", "mm", "oo"]
            value = "cc"
#elif       IK_ENABLED && D1_ENABLED
            Array = int([0,2,4,6,8,10,12,14], kind = IKG)
            value = 2_IKG
#elif       CK_ENABLED && D1_ENABLED
            Array = cmplx([0,2,4,6,8,10,12,14], -[0,2,4,6,8,10,12,14], kind = CKG)
            value = cmplx(2, 2, kind = CKG)
#elif       RK_ENABLED && D1_ENABLED
            Array = real([0,2,4,6,8,10,12,14], kind = RKG)
            value = 2._RKG
#endif
            index_ref = 2_IK

            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A value that is equal to the second element of Array, has `index = 2_IK`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_ENABLED
            Array = "acegikm"
            value = "c"
#elif       SK_ENABLED && D1_ENABLED
            Array = ["aa", "cc", "ee", "gg", "ii", "kk", "mm"]
            value = "cc"
#elif       IK_ENABLED && D1_ENABLED
            Array = int([0,2,4,6,8,10,12], kind = IKG)
            value = 2_IKG
#elif       CK_ENABLED && D1_ENABLED
            Array = cmplx([0,2,4,6,8,10,12], -[0,2,4,6,8,10,12], kind = CKG)
            value = cmplx(2, -2, kind = CKG)
#elif       RK_ENABLED && D1_ENABLED
            Array = real([0,2,4,6,8,10,12], kind = RKG)
            value = 2._RKG
#endif
            index_ref = 2_IK

            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A value that is equal to the second element of Array, has `index = 2_IK`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_ENABLED
            Array = "acegikm"
            value = "k"
#elif       SK_ENABLED && D1_ENABLED
            Array = ["aa", "cc", "ee", "gg", "ii", "kk", "mm"]
            value = "kk"
#elif       IK_ENABLED && D1_ENABLED
            Array = int([0,2,4,6,8,10,12], kind = IKG)
            value = 10_IKG
#elif       CK_ENABLED && D1_ENABLED
            Array = cmplx([0,2,4,6,8,10,12], -[0,2,4,6,8,10,12], kind = CKG)
            value = cmplx(10, -2, kind = CKG)
#elif       RK_ENABLED && D1_ENABLED
            Array = real([0,2,4,6,8,10,12], kind = RKG)
            value = 10._RKG
#endif
            index_ref = 6_IK

            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A value that is equal to the second element of Array, has `index = 2_IK`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_ENABLED
            Array = "acegikmo"
            value = "o"
#elif       SK_ENABLED && D1_ENABLED
            Array = ["aa", "cc", "ee", "gg", "ii", "kk", "mm", "oo"]
            value = "oo"
#elif       IK_ENABLED && D1_ENABLED
            Array = int([0,2,4,6,8,10,12,14], kind = IKG)
            value = 14_IKG
#elif       CK_ENABLED && D1_ENABLED
            Array = cmplx([0,2,4,6,8,10,12,14], -[0,2,4,6,8,10,12,14], kind = CKG)
            value = cmplx(14, 3, kind = CKG)
#elif       RK_ENABLED && D1_ENABLED
            Array = real([0,2,4,6,8,10,12,14], kind = RKG)
            value = 14._RKG
#endif
            index_ref = GET_SIZE(Array, kind = IK)

            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A value that is equal to the last element of Array, has `index = size(Array)`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_ENABLED
            Array = "acegikmo"
            value = "z"
#elif       SK_ENABLED && D1_ENABLED
            Array = ["aa", "cc", "ee", "gg", "ii", "kk", "mm", "oo"]
            value = "oz"
#elif       IK_ENABLED && D1_ENABLED
            Array = int([0,2,4,6,8,10,12,14], kind = IKG)
            value = 20_IKG
#elif       CK_ENABLED && D1_ENABLED
            Array = cmplx([0,2,4,6,8,10,12,14], [0,2,4,6,8,10,12,14]**2, kind = CKG)
            value = cmplx(20, -5, kind = CKG)
#elif       RK_ENABLED && D1_ENABLED
            Array = real([0,2,4,6,8,10,12,14], kind = RKG)
            value = 20._RKG
#endif
            index_ref = GET_SIZE(Array, kind = IK)

            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A value that is larger than element of Array, has `index = size(Array)`.", int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            call reset()

            Array = "acegikmo"

            value = "cg"
            index_ref = 2_IK
            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A character value """//getStr(value)//""" has `index = 2_IK`.", int(__LINE__, IK))

            value = "ca"
            index_ref = 1_IK
            call report(isLess)
            call test%assert(assertion, PROCEDURE_NAME//SK_": A character value """//getStr(value)//""" has `index = 1_IK`.", int(__LINE__, IK))
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(isLess)
            use pm_io, only: display_type
            logical(LK) , external  , optional              :: isLess
            type(display_type)                              :: disp

            if (present(isLess)) then
                index = getBin(Array, value, isLess)
            else
                index = getBin(Array, value)
            end if

            ! Report test results if needed.

            disp = display_type()
!write(*,*) getBinEnabled, present(instance), present(sorted), present(positive)
!write(*,*) Array
!write(*,*) value
!write(*,*) index
!write(*,*) index_ref
            assertion = assertion .and. index == index_ref

            if (test%traceable .and. .not. assertion) then

                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")

                call disp%show("index")
                call disp%show( index )
                call disp%show("index_ref")
                call disp%show( index_ref )
                write(test%disp%unit,"(*(g0,:,', '))") "present(isLess)    ", present(isLess)

                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP

            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GET_SIZE
#undef IS_EQUAL