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
!>  This include file contains the implementations of the tests of procedures in [pm_matrixDet](@ref pm_matrixDet).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        implicit none

        real(RK)    , parameter     :: TOLERANCE = epsilon(0._RK) * 100_IK
        real(RK)    , allocatable   :: Matrix(:,:)
        real(RK)                    :: result_ref
        real(RK)                    :: result
        real(RK)                    :: diff
        real(RK)                    :: det
#if     test_setDet_ENABLED || test_setMatDetSqrtLog_ENABLED || test_setMatDetSqrt_ENABLED
        logical(LK)                 :: failed
#endif

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Matrix = reshape(   [ 1._RK, 0._RK, 2._RK &
                            , 0._RK, 4._RK, 0._RK &
                            , 2._RK, 0._RK, 8._RK ], shape = [3,3])
        det = +16._RK
#if     test_getDet_ENABLED
        result = getMatDet(Matrix)
        result_ref = det
#elif   test_setDet_ENABLED
        call setMatDet(Matrix, result, failed)
        result_ref = det
#elif   test_getMatDetSqrtLog_ENABLED
        result = getMatDetSqrtLog(Matrix)
        result_ref = 0.5_RK * log(det)
#elif   test_setMatDetSqrtLog_ENABLED
        call setMatDetSqrtLog(Matrix, result, failed)
        result_ref = 0.5_RK * log(det)
#elif   test_getMatDetSqrt_ENABLED
        result = getMatDetSqrt(Matrix)
        result_ref = sqrt(det)
#elif   test_setMatDetSqrt_ENABLED
        call setMatDetSqrt(Matrix, result, failed)
        result_ref = sqrt(det)
#endif

#if     test_setDet_ENABLED || test_setMatDetSqrtLog_ENABLED || test_setMatDetSqrt_ENABLED
        assertion = assertion .and. .not. failed
        call test%assert(assertion, SK_"The result must be computed without failure.")
#endif

        call report()
        call test%assert(assertion, SK_"The result must be computed correctly and accurately.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_getMatDetSqrt_ENABLED || test_setMatDetSqrt_ENABLED || test_getMatDetSqrtLog_ENABLED || test_setMatDetSqrtLog_ENABLED
        Matrix = reshape(   [ 1._RK, 1._RK, 2._RK &
                            , 0._RK, 4._RK, 1._RK &
                            , 0._RK, 0._RK, 8._RK ], shape  = [3,3], order = [2,1])
        block
            use pm_matrixSymCopy, only: getMatSymFromMatUpp ! LCOV_EXCL_LINE
            det = getMatDet(getMatSymFromMatUpp(Matrix))
        end block
#if     test_getMatDetSqrtLog_ENABLED
        result = getMatDetSqrtLog(Matrix)
        result_ref = 0.5_RK * log(det)
#elif   test_setMatDetSqrtLog_ENABLED
        call setMatDetSqrtLog(Matrix, result, failed)
        result_ref = 0.5_RK * log(det)
#elif   test_getMatDetSqrt_ENABLED
        result = getMatDetSqrt(Matrix)
        result_ref = sqrt(det)
#elif   test_setMatDetSqrt_ENABLED
        call setMatDetSqrt(Matrix, result, failed)
        result_ref = sqrt(det)
#endif

#if     test_setMatDetSqrtLog_ENABLED || test_setMatDetSqrt_ENABLED
        assertion = assertion .and. .not. failed
        call test%assert(assertion, SK_"The result must be computed without failure.")
#endif

        call report()
        call test%assert(assertion, SK_"The result must be computed correctly and accurately.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_getMatDetSqrt_ENABLED || test_setMatDetSqrt_ENABLED || test_getMatDetSqrtLog_ENABLED || test_setMatDetSqrtLog_ENABLED
        Matrix = reshape(   [ 1._RK, 0._RK, 1._RK &
                            , 0._RK, 2._RK, 0._RK &
                            , 1._RK, 0._RK, 3._RK ], shape  = [3,3], order = [2,1])
        block
            use pm_matrixSymCopy, only: getMatSymFromMatUpp ! LCOV_EXCL_LINE
            det = getMatDet(getMatSymFromMatUpp(Matrix))
            if (abs(det - 4._RK) > TOLERANCE) error stop ! LCOV_EXCL_LINE
        end block
#if     test_getMatDetSqrtLog_ENABLED
        result = getMatDetSqrtLog(Matrix)
        result_ref = 0.5_RK * log(det)
#elif   test_setMatDetSqrtLog_ENABLED
        call setMatDetSqrtLog(Matrix, result, failed)
        result_ref = 0.5_RK * log(det)
#elif   test_getMatDetSqrt_ENABLED
        result = getMatDetSqrt(Matrix)
        result_ref = sqrt(det)
#elif   test_setMatDetSqrt_ENABLED
        call setMatDetSqrt(Matrix, result, failed)
        result_ref = sqrt(det)
#endif

#if     test_setMatDetSqrtLog_ENABLED || test_setMatDetSqrt_ENABLED
        assertion = assertion .and. .not. failed
        call test%assert(assertion, SK_"The result must be computed without failure.")
#endif

        call report()
        call test%assert(assertion, SK_"The result must be computed correctly and accurately.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_getDet_ENABLED || test_setDet_ENABLED
        det = -155._RK
        Matrix = reshape(   [ 10._RK, -7._RK, 0._RK &
                            , -3._RK, +2._RK, 6._RK &
                            , +5._RK, -1._RK, 5._RK ], shape  = [3,3], order = [2,1])
        result_ref = det
#if     test_getDet_ENABLED
        result = getMatDet(Matrix)
#elif   test_setDet_ENABLED
        call setMatDet(Matrix, result, failed)
        assertion = assertion .and. .not. failed
        call test%assert(assertion, SK_"The result must be computed without failure.")
#endif
        result_ref = det

        call report()
        call test%assert(assertion, SK_"The result must be computed correctly and accurately.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_getDet_ENABLED || test_setDet_ENABLED
        det = -4._RK
        Matrix = reshape(   [ -1._RK, -0._RK, -1._RK &
                            , -0._RK, -2._RK, -0._RK &
                            , -1._RK, -0._RK, -3._RK ], shape  = [3,3], order = [2,1])
        result_ref = det
#if     test_getDet_ENABLED
        result = getMatDet(Matrix)
#elif   test_setDet_ENABLED
        call setMatDet(Matrix, result, failed)
        assertion = assertion .and. .not. failed
        call test%assert(assertion, SK_"The result must be computed without failure.")
#endif
        result_ref = det

        call report()
        call test%assert(assertion, SK_"The result must be computed correctly and accurately.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_getDet_ENABLED || test_setDet_ENABLED
        det = 5070000._RK
        Matrix = reshape(   [ 17._RK, 24._RK, +1._RK, +8._RK, 15._RK &
                            , 23._RK, +5._RK, +7._RK, 14._RK, 16._RK &
                            , +4._RK, +6._RK, 13._RK, 20._RK, 22._RK &
                            , 10._RK, 12._RK, 19._RK, 21._RK, +3._RK &
                            , 11._RK, 18._RK, 25._RK, +2._RK, +9._RK ], shape  = [5,5], order = [2,1])
         result_ref = det
#if     test_getDet_ENABLED
        result = getMatDet(Matrix)
#elif   test_setDet_ENABLED
        call setMatDet(Matrix, result, failed)
        assertion = assertion .and. .not. failed
        call test%assert(assertion, SK_"The result must be computed without failure.")
#endif
        result_ref = det

        call report()
        call test%assert(assertion, SK_"The result must be computed correctly and accurately.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_setDet_ENABLED || test_setMatDetSqrtLog_ENABLED || test_setMatDetSqrt_ENABLED
        if (allocated(Matrix)) deallocate(Matrix)
        allocate(Matrix(5,5), source = 0._RK)
#if     test_setDet_ENABLED
        call setMatDet(Matrix, result, failed)
#elif   test_setMatDetSqrtLog_ENABLED
        call setMatDetSqrtLog(Matrix, result, failed)
#elif   test_setMatDetSqrt_ENABLED
        call setMatDetSqrt(Matrix, result, failed)
#endif
        assertion = assertion .and. failed
        call test%assert(assertion, SK_"The result must be computed without failure.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = (result - result_ref) / abs(result_ref)
            assertion = diff < TOLERANCE

            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "result_ref ", result_ref
                write(test%disp%unit,"(*(g0,:,', '))") "result     ", result
                write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
