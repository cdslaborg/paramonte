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
!>  This include file contains procedure implementations of the tests of [pm_distMultiNorm](@ref pm_distMultiNorm).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK)                 :: i
        integer(IK) , parameter     :: NDIM = 2_IK
        integer(IK) , parameter     :: NVEC = 20000_IK

        real(RK)    , parameter     :: TOL = 0.1_RK
        real(RK)                    :: rand(NDIM,NVEC)
        real(RK)                    :: mean(NDIM), choLow(NDIM,NDIM)
        real(RK)                    :: Mean_ref(NDIM), ChoLow_ref(NDIM,NDIM), ChoDia_ref(NDIM)
        real(RK)    , allocatable   :: diff(:)
        logical(LK)                 :: failed

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! MVN distribution
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Mean_ref = [0._RK, 0._RK]
        ChoLow_ref = reshape([1._RK, 0._RK, 0._RK, 1._RK], shape = [NDIM,NDIM])
        call setChoLow(ChoLow_ref, ChoDia_ref, failed)
        assertion = assertion .and. .not. failed
        call test%assert(assertion, desc = "The Cholesky factorization must not fail: "//getStr(ChoLow_ref))

        do i = 1_IK, NVEC
#if         getMultiNormRand_ENABLED
            rand(:,i) = getMultiNormRand(Mean_ref)
#elif       setMultiNormRand_ENABLED
            call setMultiNormRand(rand(:,i))
#endif
        end do

        call test%assert(assertion, desc = "The procedure must be able to output vector harvest.")
        call report()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! MVN distribution with non-zero mean
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Mean_ref = [-5._RK, -5._RK]
        ChoLow_ref = reshape([1._RK, 0._RK, 0._RK, 1._RK], shape = [NDIM,NDIM])
        call setChoLow(ChoLow_ref, ChoDia_ref, failed)
        assertion = assertion .and. .not. failed
        call test%assert(assertion, desc = "The Cholesky factorization must not fail: "//getStr(ChoLow_ref))

        do i = 1_IK, NVEC
#if         getMultiNormRand_ENABLED
            rand(:,i) = getMultiNormRand(Mean_ref)
#elif       setMultiNormRand_ENABLED
            call setMultiNormRand(rand(:,i), Mean_ref)
#endif
        end do

        call test%assert(assertion, desc = "The procedure must be able to output vector harvest.")
        call report()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! MVN distribution with non-identity CovMat
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Mean_ref = [0._RK, 0._RK]
        ChoLow_ref = reshape([1.5_RK, 0.6_RK, 0.6_RK, 2._RK], shape = [NDIM,NDIM])
        call setChoLow(ChoLow_ref, ChoDia_ref, failed)
        assertion = assertion .and. .not. failed
        call test%assert(assertion, desc = "The Cholesky factorization must not fail: "//getStr(ChoLow_ref))

        do i = 1_IK, NVEC
#if         getMultiNormRand_ENABLED
            rand(:,i) = getMultiNormRand(ChoLow_ref, ChoDia_ref)
#elif       setMultiNormRand_ENABLED
            call setMultiNormRand(rand(:,i), ChoLow_ref, ChoDia_ref)
#endif
        end do

        call test%assert(assertion, desc = "The procedure must be able to output vector harvest.")
        call report()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! MVN distribution with non-zero mean and non-identity CovMat
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Mean_ref = [+25._RK, -10._RK]
        ChoLow_ref = reshape([1.5_RK, 0.6_RK, 0.6_RK, 2._RK], shape = [NDIM,NDIM])
        call setChoLow(ChoLow_ref, ChoDia_ref, failed)
        assertion = assertion .and. .not. failed
        call test%assert(assertion, desc = "The Cholesky factorization must not fail: "//getStr(ChoLow_ref))

        do i = 1_IK, NVEC
#if         getMultiNormRand_ENABLED
            rand(:,i) = getMultiNormRand(Mean_ref, ChoLow_ref, ChoDia_ref)
#elif       setMultiNormRand_ENABLED
            call setMultiNormRand(rand(:,i), Mean_ref, ChoLow_ref, ChoDia_ref)
#endif
        end do

        call test%assert(assertion, desc = "The procedure must be able to output vector harvest.")
        call report()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()

            mean = getMean(rand, dim = 2_IK)

            diff = abs(mean - Mean_ref)
            assertion = assertion .and. all(diff < TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Mean_ref   ", Mean_ref
                write(test%disp%unit,"(*(g0,:,', '))") "mean       ", mean
                write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if

            call test%assert(assertion, SK_"A standard Multivariate Normal sample must have the correct mean vector within the specified tolerance. "// & ! LCOV_EXCL_LINE
                                               SK_"A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
                                               SK_"Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
                                               )

            choLow = getCovUpp(getShifted(rand, -mean, dim = 2_IK), biased = .false._LK)
            call setMatChol(choLow, choDia, info, uppDia)
            assertion = assertion .and. info == 0_IK
            call test%assert(assertion, desc = "The Cholesky factorization must not fail.")

            diff = reshape(abs(choLow - ChoLow_ref), shape = [size(choLow)])
            assertion = assertion .and. all(diff < TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "ChoLow_ref ", ChoLow_ref
                write(test%disp%unit,"(*(g0,:,', '))") "choLow     ", choLow
                write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if

            call test%assert(assertion, desc =  SK_"The output Multivariate Normal sample must have the correct Covariance and Cholesky Factorization Matrices within the specified tolerance."// & ! LCOV_EXCL_LINE
                                                SK_"A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
                                                SK_"Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
                                                )

            diff = abs(choDia - ChoDia_ref)
            assertion = assertion .and. all(diff < TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "ChoDia_ref ", ChoDia_ref
                write(test%disp%unit,"(*(g0,:,', '))") "choDia     ", choDia
                write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if

            call test%assert(assertion, desc =  SK_"The output Multivariate Normal sample must have the correct Cholesky diagonal vector within the specified tolerance."// & ! LCOV_EXCL_LINE
                                                SK_"A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
                                                SK_"Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
                                                )

        end subroutine