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
!>  This include file contains the implementations of procedures of [test_pm_complexAbs](@ref test_pm_complexAbs).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_distUnif, only: setUnifRand
        use pm_kind, only: SK, IK, LK
        use pm_val2str, only: getStr

#if     abs_ENABLED

        integer(IK)     , parameter :: NP = 5_IK
        character(*, SK), parameter :: PROCEDURE_NAME = "abs()"
        complex(CKC)    :: Matrix(NP,NP), MatrixABS(NP,NP)
        integer(IK)     :: i, j

        assertion = .true._LK
        call setUnifRand(Matrix(1,1))
        MatrixABS(1,1) = abs(Matrix(1,1))

        assertion = assertion .and. abs(Matrix(1,1)%re) == MatrixABS(1,1)%re
        call report(1_IK, 1_IK, int(__LINE__, IK))

        assertion = assertion .and. abs(Matrix(1,1)%im) == MatrixABS(1,1)%im
        call report(1_IK, 1_IK, int(__LINE__, IK))

        call setUnifRand(Matrix)
        MatrixABS = abs(Matrix)

        assertion = assertion .and. all(abs(Matrix%re) == MatrixABS%re)
        do j = 1, NP
            do i = 1, NP
                call report(1_IK, 1_IK, int(__LINE__, IK))
            end do
        end do

        assertion = assertion .and. all(abs(Matrix%im) == MatrixABS%im)
        do j = 1, NP
            do i = 1, NP
                call report(1_IK, 1_IK, int(__LINE__, IK))
            end do
        end do

    contains

        subroutine report(i,j, line)
            integer(IK), intent(in) :: i, j, line
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "abs(Matrix(i,j)%re)", abs(Matrix(i,j)%re)
                write(test%disp%unit,"(*(g0,:,', '))") "abs(Matrix(i,j)%im)", abs(Matrix(i,j)%im)
                write(test%disp%unit,"(*(g0,:,', '))") "MatrixABS(i,j%re)  ", MatrixABS(i,j)%re
                write(test%disp%unit,"(*(g0,:,', '))") "MatrixABS(i,j%im)  ", MatrixABS(i,j)%im
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, PROCEDURE_NAME//SK_" The elemental absolute of complex values must be correctly computed.", line)
        end subroutine

#else

#error  "Unrecognized interface."

#endif

#undef COMPARE
