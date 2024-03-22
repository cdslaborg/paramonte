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
!>  This include file contains procedure implementations of the tests of [getNormLogPDF](@ref pm_distNorm::getNormLogPDF).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getNormLogPDF_ENABLED || setNormLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RK)    , parameter     :: TOL = epsilon(1._RK) * 100_IK
        real(RK)    , parameter     :: Point(*) = [-3.5_RK, -1._RK, 0._RK, .95_RK, 1.5_RK]
        real(RK)    , allocatable   :: LogPDF_ref(:)
        real(RK)    , allocatable   :: logPDF(:)
        real(RK)    , allocatable   :: diff(:)
        real(RK)    , allocatable   :: logInvStd
        real(RK)    , allocatable   :: invStd
        real(RK)    , allocatable   :: mean
        integer(IK)                 :: i

        assertion = .true._LK

        ! Compute the logPDF with mean and invStd

        mean = 2.5_RK
        invStd = 1._RK / 3.3_RK
        logInvStd = log(invStd)
        LogPDF_ref =    [ -3.76575356366057836759969238586510697_RK & ! LCOV_EXCL_LINE
                        , -2.67530360957426064491833334086969790_RK & ! LCOV_EXCL_LINE
                        , -2.39982151591034879918830579266033172_RK & ! LCOV_EXCL_LINE
                        , -2.22316862334836532811392562737107537_RK & ! LCOV_EXCL_LINE
                        , -2.15877468395442593417453168797713593_RK & ! LCOV_EXCL_LINE
                        ]

        logPDF = [( getNormLogPDF(Point(i), mean, invStd, logInvStd), i = 1, size(Point) )]
        call report()
        call test%assert(assertion, desc = "The logPDF must be computed correctly with an input mean and invStd.")

        logPDF = getNormLogPDF(Point-mean, invStd, logInvStd)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input invStd.")

        ! Compute the logPDF with mean and invStd

        mean = 2.5_RK
        invStd = 1._RK / 3.3_RK
        logInvStd = log(invStd)
        LogPDF_ref = getNormLogPDF(Point, mean, invStd, logInvStd)

        logPDF = [( getNormLogPDF(Point(i) - mean, invStd, logInvStd), i = 1, size(Point) )]
        call report()
        call test%assert(assertion, desc = "The logPDF must be computed correctly with an input invStd.")

        logPDF = getNormLogPDF(Point, mean, invStd, logInvStd)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input mean and invStd.")

        ! Compute the logPDF with mean

        mean = 2.5_RK
        invStd = 1._RK
        logInvStd = 0._RK
        LogPDF_ref = getNormLogPDF(Point, mean, invStd, logInvStd)

        logPDF = [( getNormLogPDF(Point(i), mean), i = 1, size(Point) )]
        call report()
        call test%assert(assertion, desc = "The logPDF must be computed correctly with an input mean.")

        logPDF = getNormLogPDF(Point, mean)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input mean.")

        ! Compute the logPDF with mean

        mean = 0._RK
        invStd = 1._RK
        LogPDF_ref = getNormLogPDF(Point, mean, invStd, logInvStd)

        logPDF = [( getNormLogPDF(Point(i)), i = 1, size(Point) )]
        call report()
        call test%assert(assertion, desc = "The logPDF must be computed correctly.")

        logPDF = getNormLogPDF(Point)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs.")

    contains

        subroutine report()
            diff = abs(logPDF - LogPDF_ref)
            assertion = all(diff < TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "LogPDF_ref ", LogPDF_ref
                write(test%disp%unit,"(*(g0,:,', '))") "logPDF     ", logPDF
                write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getNormCDF_ENABLED || setNormCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RK)    , parameter     :: TOL = epsilon(1._RK) * 100_IK
        real(RK)    , parameter     :: Point(*) = [-3.5_RK, -1._RK, 0._RK, .95_RK, 1.5_RK]
        real(RK)    , allocatable   :: CDF_ref(:)
        real(RK)    , allocatable   :: CDF(:)
        real(RK)    , allocatable   :: diff(:)
        real(RK)    , allocatable   :: invStd
        real(RK)    , allocatable   :: mean
        integer(IK)                 :: i

        assertion = .true._LK

        ! Compute the CDF with mean and invStd

        mean = 2.5_RK
        invStd = 1._RK / 3.3_RK
        CDF_ref =   [ 0.345181739972076390427473231450096989E-1_RK  & ! LCOV_EXCL_LINE
                    , 0.144434483561982577092518129626705819_RK     & ! LCOV_EXCL_LINE
                    , 0.224352498756817662874749900906318120_RK     & ! LCOV_EXCL_LINE
                    , 0.319285766588537302588523434658987546_RK     & ! LCOV_EXCL_LINE
                    , 0.380933384094123902018957824101812120_RK     & ! LCOV_EXCL_LINE
                    ]

        CDF = [( getNormCDF(Point(i), mean, invStd), i = 1, size(Point) )]
        call report()
        call test%assert(assertion, desc = "The CDF must be computed correctly with an input mean and invStd.")

        CDF = getNormCDF(Point, mean, invStd)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input mean and invStd.")

        CDF = getNormCDF(Point-mean, invStd = invStd)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input invStd.")

        CDF = getNormCDF(Point*invStd, mean = mean*invStd)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input mean.")

        CDF = getNormCDF((Point-mean)*invStd)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with no input parameters.")


        CDF = getNormCDF(Point, spread(mean,1,size(Point)), spread(invStd,1,size(Point)))
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input mean and invStd vectors.")

        CDF = getNormCDF(Point-mean, invStd = spread(invStd,1,size(Point)))
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input invStd vector.")

        CDF = getNormCDF(Point*invStd, mean = mean*spread(invStd,1,size(Point)))
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input mean vector.")

        ! Compute the CDF with mean

        mean = 2.5_RK
        invStd = 1._RK
        CDF_ref = getNormCDF(Point, mean, invStd)

        CDF = [( getNormCDF(Point(i), mean), i = 1, size(Point) )]
        call report()
        call test%assert(assertion, desc = "The CDF must be computed correctly with an input mean.")

        CDF = getNormCDF(Point, mean)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with an input mean.")

        ! Compute the CDF with mean

        mean = 0._RK
        invStd = 1._RK
        CDF_ref = getNormCDF(Point, mean, invStd)

        CDF = [( getNormCDF(Point(i)), i = 1, size(Point) )]
        call report()
        call test%assert(assertion, desc = "The CDF must be computed correctly.")

        CDF = getNormCDF(Point)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs.")

        CDF = getNormCDF(Point, mean)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with mean present.")

        CDF = getNormCDF(Point, invStd = invStd)
        call report()
        call test%assert(assertion, desc = "The procedure must be able to act on 1D array inputs with invStd present.")

    contains

        subroutine report()
            diff = abs(CDF - CDF_ref)
            assertion = all(diff < TOL)

            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "CDF_ref    ", CDF_ref
                write(test%disp%unit,"(*(g0,:,', '))") "CDF        ", CDF
                write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getNormQuan_ENABLED || setNormQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getNormRand_ENABLED || setNormRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getNormRand_ENABLED
        integer(IK)             :: i
#endif
        integer(IK) , parameter :: NP1 = 40_IK
        integer(IK) , parameter :: NP2 = 50_IK
        integer(IK) , parameter :: NP3 = 60_IK
        real(RK)    , parameter :: TOL = 0.05_RK
        real(RK)                :: std_ref
        real(RK)                :: mean_ref
        real(RK)                :: Sample(NP1,NP2,NP3)
        real(RK)                :: mean, std, diff
        assertion = .true._LK

        ! Standard normal distribution

        mean_ref = 0._RK
        std_ref = 1._RK

#if     getNormRand_ENABLED
        Sample(1,1,1) = getNormRand(0._RK)
#elif   setNormRand_ENABLED
        call setNormRand(Sample(1,1,1))
#endif
        call test%assert(assertion, desc = "The procedure must be able to output scalar harvest.")

#if     getNormRand_ENABLED
        Sample(:,1,1) = getNormRand([(0._RK, i = 1, NP1)])
#elif   setNormRand_ENABLED
        call setNormRand(Sample(:,1,1))
#endif
        call test%assert(assertion, desc = "The procedure must be able to output 1D array harvest.")

#if     getNormRand_ENABLED
        Sample(:,:,1) = getNormRand(reshape([(0._RK, i = 1, NP1*NP2)], shape = [NP1,NP2]))
#elif   setNormRand_ENABLED
        call setNormRand(Sample(:,:,1))
#endif
        call test%assert(assertion, desc = "The procedure must be able to output 2D array harvest.")

#if     getNormRand_ENABLED
        Sample(:,:,:) = getNormRand(reshape([(0._RK, i = 1, NP1*NP2*NP3)], shape = [NP1,NP2,NP3]))
#elif   setNormRand_ENABLED
        call setNormRand(Sample(:,:,:))
#endif
        call test%assert(assertion, desc = "The procedure must be able to output 3D array harvest.")

        mean = getMean(reshape(Sample,[size(Sample)]))
        diff = abs(mean - mean_ref)
        assertion = diff < TOL

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "mean_ref   ", mean_ref
            write(test%disp%unit,"(*(g0,:,', '))") "mean       ", mean
            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
            write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        call test%assert(assertion, desc = "A standard Normal sample must have a mean of `0._RK`. "// & ! LCOV_EXCL_LINE
                                           "A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
                                           "Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
                                           )

        std = sqrt( getVar( reshape(Sample,[size(Sample)]), mean = mean, biased = .false._LK ) )
        diff = abs(std - std_ref)
        assertion = diff < TOL

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "mean_ref   ", mean_ref
            write(test%disp%unit,"(*(g0,:,', '))") "mean       ", mean
            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
            write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        call test%assert(assertion, desc = "A standard Normal sample must have a standard deviation of `1._RK`. "// & ! LCOV_EXCL_LINE
                                           "A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
                                           "Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
                                           )

        ! normal distribution with mean

        mean_ref = -2.33_RK
        std_ref = 1._RK

#if     getNormRand_ENABLED
        Sample(1,1,1) = getNormRand(mean = mean_ref)
#elif   setNormRand_ENABLED
        call setNormRand(Sample(1,1,1), mean = mean_ref)
#endif
        call test%assert(assertion, desc = "The procedure must be able to output scalar harvest with an input mean.")

#if     getNormRand_ENABLED
        Sample(:,1,1) = getNormRand(mean = [(mean_ref, i = 1, NP1)])
#elif   setNormRand_ENABLED
        call setNormRand(Sample(:,1,1), mean = mean_ref)
#endif
        call test%assert(assertion, desc = "The procedure must be able to output 1D array harvest with an input mean.")

#if     getNormRand_ENABLED
        Sample(:,:,1) = getNormRand(mean = reshape([(mean_ref, i = 1, NP1*NP2)], shape = [NP1, NP2]))
#elif   setNormRand_ENABLED
        call setNormRand(Sample(:,:,1), mean = mean_ref)
#endif
        call test%assert(assertion, desc = "The procedure must be able to output 2D array harvest with an input mean.")

#if     getNormRand_ENABLED
        Sample(:,:,:) = getNormRand(mean = reshape([(mean_ref, i = 1, NP1*NP2*NP3)], shape = [NP1, NP2, NP3]))
#elif   setNormRand_ENABLED
        call setNormRand(Sample(:,:,:), mean = mean_ref)
#endif
        call test%assert(assertion, desc = "The procedure must be able to output 3D array harvest with an input mean.")

        mean = getMean(reshape(Sample,[size(Sample)]))
        diff = abs(mean - mean_ref)
        assertion = diff < TOL

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "mean_ref   ", mean_ref
            write(test%disp%unit,"(*(g0,:,', '))") "mean       ", mean
            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
            write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        call test%assert(assertion, desc = "A Normal sample with a given input mean must have the same mean. "// & ! LCOV_EXCL_LINE
                                           "A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
                                           "Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
                                           )

        std = sqrt( getVar( reshape(Sample,[size(Sample)]), mean = mean, biased = .false._LK ) )
        diff = abs(std - std_ref)
        assertion = diff < TOL

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "mean_ref   ", mean_ref
            write(test%disp%unit,"(*(g0,:,', '))") "mean       ", mean
            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
            write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        call test%assert(assertion, desc = "A Normal sample with missing input std must have a std of `1._RK`. "// & ! LCOV_EXCL_LINE
                                           "A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
                                           "Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
                                           )

        ! normal distribution with std and without mean

!        mean_ref = 0._RK
!        std_ref = 2._RK
!
!        call setNormRand(Sample(1,1,1), std = std_ref)
!        call test%assert(assertion, desc = "The procedure must be able to output scalar harvest with an input std.")
!
!        call setNormRand(Sample(:,1,1), std = std_ref)
!        call test%assert(assertion, desc = "The procedure must be able to output 1D array harvest with an input std.")
!
!        call setNormRand(Sample(:,:,1), std = std_ref)
!        call test%assert(assertion, desc = "The procedure must be able to output 2D array harvest with an input std.")
!
!        call setNormRand(Sample(:,:,:), std = std_ref)
!        call test%assert(assertion, desc = "The procedure must be able to output 3D array harvest with an input std.")
!
!        mean = getMean(reshape(Sample,[size(Sample)]))
!        diff = abs(mean - mean_ref)
!        assertion = diff < TOL
!
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "mean_ref   ", mean_ref
!            write(test%disp%unit,"(*(g0,:,', '))") "mean       ", mean
!            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
!            write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
!            write(test%disp%unit,"(*(g0,:,', '))")
!            ! LCOV_EXCL_STOP
!        end if
!
!        call test%assert(assertion, desc = "A Normal sample with a given input std must have the same std with zero mean. "// & ! LCOV_EXCL_LINE
!                                           "A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
!                                           "Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
!                                           )
!
!        std = sqrt( getVar( reshape(Sample,[size(Sample)]), mean = mean, biased = .false._LK ) )
!        diff = abs(std - std_ref)
!        assertion = diff < TOL
!
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0,:,', '))")
!            write(test%disp%unit,"(*(g0,:,', '))") "mean_ref   ", mean_ref
!            write(test%disp%unit,"(*(g0,:,', '))") "mean       ", mean
!            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
!            write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
!            write(test%disp%unit,"(*(g0,:,', '))")
!            ! LCOV_EXCL_STOP
!        end if
!
!        call test%assert(assertion, desc = "A Normal sample with missing input mean must have a mean of `0._RK`. "// & ! LCOV_EXCL_LINE
!                                           "A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
!                                           "Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
!                                           )

        ! normal distribution with mean and std

        mean_ref = 5._RK
        std_ref = 3.25_RK

#if     getNormRand_ENABLED
        Sample(1,1,1) = getNormRand(mean = mean_ref, std = std_ref)
#elif   setNormRand_ENABLED
        call setNormRand(Sample(1,1,1), mean = mean_ref, std = std_ref)
#endif
        call test%assert(assertion, desc = "The procedure must be able to output scalar harvest with an input mean and std.")

#if     getNormRand_ENABLED
        Sample(:,1,1) = getNormRand(mean = [(mean_ref, i = 1, NP1)], std = std_ref)
#elif   setNormRand_ENABLED
        call setNormRand(Sample(:,1,1), mean = mean_ref, std = std_ref)
#endif
        call test%assert(assertion, desc = "The procedure must be able to output 1D array harvest with an input mean and std.")

#if     getNormRand_ENABLED
        Sample(:,:,1) = getNormRand(mean = reshape([(mean_ref, i = 1, NP1*NP2)], shape = [NP1, NP2]), std = std_ref)
#elif   setNormRand_ENABLED
        call setNormRand(Sample(:,:,1), mean = mean_ref, std = std_ref)
#endif
        call test%assert(assertion, desc = "The procedure must be able to output 2D array harvest with an input mean and std.")

#if     getNormRand_ENABLED
        Sample(:,:,:) = getNormRand(mean = reshape([(mean_ref, i = 1, NP1*NP2*NP3)], shape = [NP1, NP2, NP3]), std = std_ref)
#elif   setNormRand_ENABLED
        call setNormRand(Sample(:,:,:), mean = mean_ref, std = std_ref)
#endif
        call test%assert(assertion, desc = "The procedure must be able to output 3D array harvest with an input mean and std.")

        mean = getMean(reshape(Sample,[size(Sample)]))
        diff = abs(mean - mean_ref)
        assertion = diff < TOL

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "mean_ref   ", mean_ref
            write(test%disp%unit,"(*(g0,:,', '))") "mean       ", mean
            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
            write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        call test%assert(assertion, desc = "A Normal sample with a given input mean must have the same mean. "// & ! LCOV_EXCL_LINE
                                           "A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
                                           "Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
                                           )

        std = sqrt( getVar( reshape(Sample,[size(Sample)]), mean = mean, biased = .false._LK ) )
        diff = abs(std - std_ref)
        assertion = diff < TOL

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "mean_ref   ", mean_ref
            write(test%disp%unit,"(*(g0,:,', '))") "mean       ", mean
            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
            write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        call test%assert(assertion, desc = "A Normal sample with a given input std must have the same std. "// & ! LCOV_EXCL_LINE
                                           "A failure of this test could be due to the stochastic nature of random number generation. "// & ! LCOV_EXCL_LINE
                                           "Repeating the test may resolve the failure. If the problem persists, increase the `TOL` parameter in the test. " & ! LCOV_EXCL_LINE
                                           )

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getNormEntropy_ENABLED || setNormEntropy_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RK)    , parameter     :: TOL = epsilon(1._RK) * 100_IK
        real(RK)    , parameter     :: Variance(*) = [0.1_RK, 0.5_RK, 10.8_RK]
        real(RK)    , allocatable   :: Entropy_ref(:)
        real(RK)    , allocatable   :: Entropy(:)
        real(RK)    , allocatable   :: diff(:)
        integer(IK)                 :: i

        assertion = .true._LK

        Entropy_ref = 0.5_RK * log(2*acos(-1._RK)*Variance) + 0.5_RK
        Entropy = [( getNormEntropy(logVar = log(Variance(i))), i = 1, size(Variance) )]
        diff = abs(Entropy - Entropy_ref)
        assertion = all(diff < TOL)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "Entropy_ref    ", Entropy_ref
            write(test%disp%unit,"(*(g0,:,', '))") "Entropy        ", Entropy
            write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
            write(test%disp%unit,"(*(g0,:,', '))") "TOL            ", TOL
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        call test%assert(assertion, desc = "The entropy must be computed correctly.")

        Entropy = getNormEntropy(logVar = log(Variance))
        diff = abs(Entropy - Entropy_ref)
        assertion = all(diff < TOL)

        call test%assert(assertion, desc = "The procedure must be able to act on array inputs.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getNormFisher_ENABLED || setNormFisher_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RK)    , parameter     :: TOL = epsilon(1._RK) * 100_IK
        real(RK)    , parameter     :: Point(*) = [-3.5_RK, -1._RK, 0._RK, .95_RK, 1.5_RK]
        real(RK)    , allocatable   :: Fisher_ref(:,:)
        real(RK)    , allocatable   :: Fisher(:,:)
        real(RK)    , allocatable   :: diff(:,:)
        real(RK)    , allocatable   :: varInv

        assertion = .true._LK

        ! Compute the Fisher with mean and varInv

        varInv = 1._RK / 3.3_RK**2
        Fisher_ref = reshape(   [ varInv & ! LCOV_EXCL_LINE
                                , 0._RK & ! LCOV_EXCL_LINE
                                , 0._RK & ! LCOV_EXCL_LINE
                                , 2 * varInv & ! LCOV_EXCL_LINE
                                ], shape = [2,2])
        Fisher = getNormFisher(varInv)
        diff = abs(Fisher - Fisher_ref)
        assertion = all(diff < TOL)

        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "Fisher_ref ", Fisher_ref
            write(test%disp%unit,"(*(g0,:,', '))") "Fisher     ", Fisher
            write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
            write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
            write(test%disp%unit,"(*(g0,:,', '))")
            ! LCOV_EXCL_STOP
        end if

        call test%assert(assertion, desc = "The Fisher must be computed correctly with an input varInv.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getNormKLD_ENABLED || setNormKLD_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) , parameter     :: NP = 5_IK
        real(RK)    , parameter     :: TOL = epsilon(1._RK) * 100_IK
        real(RK)    , parameter     :: MeanP(NP) = [0._RK, .1_RK, 1.0_RK, -1.0_RK,  1.5_RK]
        real(RK)    , parameter     :: VarP(NP)  = [2._RK, 1._RK, 2.0_RK, +1.5_RK,  10._RK]
        real(RK)    , parameter     :: MeanQ(NP) = [0._RK, .3_RK, 1.5_RK, -1.3_RK,  2.0_RK]
        real(RK)    , parameter     :: VarQ(NP)  = [3._RK, 1._RK, 3.0_RK, +2.5_RK, 11.5_RK]
        real(RK)    , parameter     :: MeanDiffSq(NP) = (MeanP - MeanQ)**2
        real(RK)    , allocatable   :: KLD_ref(:)
        real(RK)    , allocatable   :: diff(:)
        real(RK)    , allocatable   :: KLD(:)
        integer(IK)                 :: i

        assertion = .true._LK

        ! Compute the KLD with MeanDiffSq and variances

        KLD_ref =   [ 0.360658873874155243223398910655078833E-1_RK &
                    , 0.199999999999999999999999999999999974E-1_RK &
                    , 0.777325540540821909890065577321745419E-1_RK &
                    , 0.734128118829953416027570481518310383E-1_RK &
                    , 0.155331451006228269466341930512189660E-1_RK ]
        KLD = [( getNormKLD(MeanDiffSq(i), VarP(i), VarQ(i)), i = 1, NP )]
        call report()
        call test%assert(assertion, desc = "The KLD must be computed correctly with input scalar MeanDiffSq, VarP, and VarQ.")

        KLD = getNormKLD(MeanDiffSq, VarP, VarQ)
        call report()
        call test%assert(assertion, desc = "The KLD must be computed correctly with input vector MeanDiffSq, VarP, and VarQ.")

        ! Compute the KLD with MeanDiffSq

        KLD_ref = getNormKLD(MeanDiffSq, 1._RK, 1._RK)

        KLD = [( getNormKLD(MeanDiffSq(i)), i = 1, NP )]
        call report()
        call test%assert(assertion, desc = "The KLD must be computed correctly with input scalar MeanDiffSq.")

        KLD = getNormKLD(MeanDiffSq)
        call report()
        call test%assert(assertion, desc = "The KLD must be computed correctly with input vector MeanDiffSq.")

        ! Compute the KLD with variances

        KLD_ref = getNormKLD(0._RK, VarP, VarQ)

        KLD = [( getNormKLD(VarP(i), VarQ(i)), i = 1, NP )]
        call report()
        call test%assert(assertion, desc = "The KLD must be computed correctly with input scalar VarP, and VarQ.")

        KLD = getNormKLD(VarP, VarQ)
        call report()
        call test%assert(assertion, desc = "The KLD must be computed correctly with input vector VarP, and VarQ.")

    contains

        subroutine report()
            diff = abs(KLD - KLD_ref)
            assertion = all(diff < TOL)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "KLD_ref    ", KLD_ref
                write(test%disp%unit,"(*(g0,:,', '))") "KLD        ", KLD
                write(test%disp%unit,"(*(g0,:,', '))") "diff       ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL        ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif