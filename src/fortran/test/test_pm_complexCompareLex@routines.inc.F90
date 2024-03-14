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
!>  This include file contains the implementations of the tests of procedures of [test_pm_complexCompareLex](@ref test_pm_complexCompareLex).
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_distUnif, only: setUnifRand
        use pm_kind, only: SK, IK, LK
        use pm_val2str, only: getStr

#if     islexless_CK_ENABLED
#define COMPARES_WITH <
        character(*, SK), parameter :: PROCEDURE_NAME = "islexless()"
#elif   islexleq_CK_ENABLED
#define COMPARES_WITH <=
        character(*, SK), parameter :: PROCEDURE_NAME = "islexleq()"
#elif   islexmeq_CK_ENABLED
#define COMPARES_WITH >=
        character(*, SK), parameter :: PROCEDURE_NAME = "islexmeq()"
#elif   islexmore_CK_ENABLED
#define COMPARES_WITH >
        character(*, SK), parameter :: PROCEDURE_NAME = "islexmore()"
#else
#error  "Unrecognized interface."
#endif

        integer(IK)     , parameter :: NP = 5_IK
        integer(IK)     :: i, j
        logical(LK)     :: result(NP,NP), result_def(NP,NP)
        complex(CKC)    :: mat1(NP,NP)
        complex(CKC)    :: mat2(NP,NP)

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        do j = 1, NP
            do i = 1, NP
                call setUnifRand(mat1(i,j))
                call setUnifRand(mat2(i,j))
                result(i,j) = mat1(i,j) COMPARES_WITH mat2(i,j)
                result_def(i,j) = iscomparable(mat1(i,j), mat2(i,j))
                assertion = assertion .and. (result(i,j) .eqv. result_def(i,j))
            end do
        end do
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must return correctly compare scalar complex values.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        do j = 1, NP
            call setUnifRand(mat1(:,j))
            call setUnifRand(mat2(:,j))
            result(:,j) = mat1(:,j) COMPARES_WITH mat2(:,j)
            result_def(:,j) = iscomparable(mat1(:,j), mat2(:,j))
            assertion = assertion .and. all(result(:,j) .eqv. result_def(:,j))
        end do
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must return correctly compare vector complex values.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setUnifRand(mat1(:,:))
        call setUnifRand(mat2(:,:))
        result(:,:) = mat1(:,:) COMPARES_WITH mat2(:,:)
        result_def(:,:) = iscomparable(mat1(:,:), mat2(:,:))
        assertion = assertion .and. all(result(:,:) .eqv. result_def(:,:))
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must return correctly compare matrix complex values.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure elemental function iscomparable(val1, val2) result(comparable)
            complex(CKC), intent(in) :: val1, val2
            logical(LK) :: comparable
#if         islexless_CK_ENABLED
            comparable = (val1%re < val2%re) .or. (val1%re == val2%re .and. val1%re < val2%re)
#elif       islexleq_CK_ENABLED
            comparable = (val1%re < val2%re) .or. (val1%re == val2%re .and. val1%re <= val2%re)
#elif       islexmeq_CK_ENABLED
            comparable = (val1%re > val2%re) .or. (val1%re == val2%re .and. val1%re >= val2%re)
#elif       islexmore_CK_ENABLED
            comparable = (val1%re > val2%re) .or. (val1%re == val2%re .and. val1%re > val2%re)
#else
#error      "Unrecognized interface."
#endif
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                loopOverRows: do j = 1, NP
                    do i = 1, NP
                        if (result(i,j) .neqv. result_def(i,j)) then
                            write(test%disp%unit,"(*(g0,:,', '))")
                            write(test%disp%unit,"(*(g0,:,', '))") "mat1(i,j)   ", mat1(i,j)
                            write(test%disp%unit,"(*(g0,:,', '))") "mat2(i,j)   ", mat2(i,j)
                            write(test%disp%unit,"(*(g0,:,', '))")
                            exit loopOverRows
                        end if
                    end do
                end do loopOverRows
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef COMPARES_WITH_DEF
#undef COMPARES_WITH