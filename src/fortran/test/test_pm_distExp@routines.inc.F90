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
!>  This include file contains procedure implementations of the tests of [pm_distExp](@ref pm_distExp).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getExpLogPDF_ENABLED || setExpLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK_ENABLED
        real(RKC)   , parameter :: TOL = epsilon(0._RKC) * 10
        real(RKC)   , parameter :: ZERO = 0._RKC, ONE = 1._RKC
        real(RKC)               :: x, mu, invSigma, logInvSigma
        real(RKC)               :: logPDF, logPDF_ref
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: itry
        assertion = .true._LK
        do itry = 1, 100
            mu = getChoice([ZERO, getUnifRand(-10._RKC, 10._RKC)])
            invSigma = getChoice([ONE, getUnifRand(0.1_RKC, 10._RKC)])
            logInvSigma = log(invSigma)
            x = getUnifRand(mu, mu + 10._RKC)
            logPDF_ref = getLogPDF_ref(x, mu, invSigma)
#if         getExpLogPDF_ENABLED
            logPDF = getExpLogPDF(x, mu, invSigma)
            call report(__LINE__, mu, invSigma)
            if (mu == ZERO) then
                logPDF = getExpLogPDF(x, invSigma = invSigma)
                call report(__LINE__, invSigma = invSigma)
            end if
            if (invSigma == ONE) then
                logPDF = getExpLogPDF(x, mu)
                call report(__LINE__, mu)
            end if
            if (mu == ZERO .and. invSigma == ONE) then
                logPDF = getExpLogPDF(x)
                call report(__LINE__)
            end if
#elif       setExpLogPDF_ENABLED
            call setExpLogPDF(logPDF, x, mu, invSigma, logInvSigma)
            call report(__LINE__, mu, invSigma)
            if (mu == ZERO) then
                call setExpLogPDF(logPDF, x, invSigma, logInvSigma)
                call report(__LINE__, invSigma = invSigma)
                if (invSigma == ONE) then
                    call setExpLogPDF(logPDF, x)
                    call report(__LINE__, mu)
                end if
            end if
#else
#error      "Unrecognized interface."
#endif
        end do

    contains

        pure elemental function getLogPDF_ref(x, mu, invSigma) result(logPDF)
            real(RKC), intent(in) :: x, mu, invSigma
            real(RKC) :: logPDF
            logPDF = logInvSigma - (x - mu) * invSigma
        end function

        subroutine report(line, mu, invSigma)
            integer(IK) , intent(in)            :: line
            real(RKC)   , intent(in), optional  :: mu, invSigma
            assertion = assertion .and. isClose(logPDF, logPDF_ref, abstol = TOL)
            call test%assert(assertion, SK_"The PDF must be computed correctly.", int(line, IK))
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("x")
                call test%disp%show( x )
                call test%disp%show("present(mu)")
                call test%disp%show( present(mu) )
                if (present(mu)) then
                call test%disp%show("mu")
                call test%disp%show( mu )
                end if
                call test%disp%show("present(invSigma)")
                call test%disp%show( present(invSigma) )
                if (present(invSigma)) then
                call test%disp%show("logInvSigma")
                call test%disp%show( logInvSigma )
                call test%disp%show("invSigma")
                call test%disp%show( invSigma )
                end if
                call test%disp%show("logPDF_ref")
                call test%disp%show( logPDF_ref )
                call test%disp%show("logPDF")
                call test%disp%show( logPDF )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getExpCDF_ENABLED || setExpCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK_ENABLED
        real(RKC)   , parameter :: TOL = epsilon(0._RKC) * 10
        real(RKC)   , parameter :: ZERO = 0._RKC, ONE = 1._RKC
        real(RKC)               :: invSigma, mu, x
        real(RKC)               :: cdf, cdf_ref
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: itry
        assertion = .true._LK
        do itry = 1, 100
            mu = getChoice([ZERO, getUnifRand(-10._RKC, 10._RKC)])
            invSigma = getChoice([ONE, getUnifRand(1._RKC, 10._RKC)])
            x = getUnifRand(mu, mu + 10._RKC)
            cdf_ref = getCDF_ref(x, mu, invSigma)
#if         getExpCDF_ENABLED
            cdf = getExpCDF(x, mu, invSigma)
            call report(__LINE__, mu, invSigma)
            if (mu == ZERO) then
                cdf = getExpCDF(x, invSigma = invSigma)
                call report(__LINE__, invSigma = invSigma)
            end if
            if (invSigma == ONE) then
                cdf = getExpCDF(x, mu)
                call report(__LINE__, mu)
            end if
            if (mu == ZERO .and. invSigma == ONE) then
                cdf = getExpCDF(x)
                call report(__LINE__)
            end if
#elif       setExpCDF_ENABLED
            call setExpCDF(cdf, x, mu, invSigma)
            call report(__LINE__, mu, invSigma)
            if (mu == ZERO) then
                call setExpCDF(cdf, x, invSigma)
                call report(__LINE__, invSigma = invSigma)
                if (invSigma == ONE) then
                    call setExpCDF(cdf, x)
                    call report(__LINE__, mu)
                end if
            end if
#else
#error      "Unrecognized interface."
#endif
        end do

    contains

        pure elemental function getCDF_ref(x, mu, invSigma) result(cdf)
            real(RKC), intent(in) :: x, mu, invSigma
            real(RKC) :: cdf
            cdf = 1._RKC - exp(-(x - mu) * invSigma)
        end function

        subroutine report(line, mu, invSigma)
            integer(IK) , intent(in)            :: line
            real(RKC)   , intent(in), optional  :: mu, invSigma
            assertion = assertion .and. isClose(cdf, cdf_ref, abstol = TOL)
            call test%assert(assertion, SK_"The CDF must be computed correctly.", int(line, IK))
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("x")
                call test%disp%show( x )
                call test%disp%show("present(mu)")
                call test%disp%show( present(mu) )
                if (present(mu)) then
                call test%disp%show("mu")
                call test%disp%show( mu )
                end if
                call test%disp%show("present(invSigma)")
                call test%disp%show( present(invSigma) )
                if (present(invSigma)) then
                call test%disp%show("invSigma")
                call test%disp%show( invSigma )
                end if
                call test%disp%show("cdf_ref")
                call test%disp%show( cdf_ref )
                call test%disp%show("cdf")
                call test%disp%show( cdf )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getExpRand_ENABLED || setExpRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_kind, only: IK
        use pm_val2str, only: getStr
        use pm_option, only: getOption
        use pm_sampleMean, only: getMean
        use pm_distUnif, only: getUnifRand
        use pm_distUnif, only: setUnifRand

        integer(IK) :: i
        integer(IK) , parameter :: NP = 10000_IK
        real(RKC)   , parameter :: TOL = 0.1_RKC
        real(RKC)   , parameter :: SIGMA_DEF = 1._RKC, MU_DEF = 0._RKC
        real(RKC)               :: rand(NP), diff, reltol

        diff = huge(0._RKC)
        reltol = -huge(0._RKC)
        assertion = .true._LK

#if     setExpRand_ENABLED
        call runTestsWith()
#endif
        do i = 1, 5
            call runTestsWith(sigma = getUnifRand(.5_RKC, 5._RKC))
            call runTestsWith(sigma = getUnifRand(.5_RKC, 5._RKC), mu = getUnifRand(-.5_RKC, 5._RKC))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(sigma, mu)
            real(RKC), intent(in), optional :: sigma, mu
            real(RKC) :: sigmac

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         getExpRand_ENABLED
            integer :: i
            if (present(mu)) then
                do i = 1, size(rand)
                    rand(i) = getExpRand(sigma, mu)
                end do
            else
                do i = 1, size(rand)
                    rand(i) = getExpRand(sigma)
                end do
            end if
#elif       setExpRand_ENABLED
            call setUnifRand(rand)
            if (present(sigma) .and. present(mu)) then
                call setExpRand(rand, sigma, mu)
            elseif (present(sigma)) then
                call setExpRand(rand, sigma)
            else
                call setExpRand(rand)
            end if
#else
#error      "Unrecognized interface."
#endif
            assertion = assertion .and. all(getOption(MU_DEF, mu) <= rand)
            call report(sigma, mu)
            call test%assert(assertion, SK_"The condition `all(mu_ref <= rand)` must hold. Repeating the test may resolve the failure. present(sigma), present(mu) = "//getStr([present(sigma), present(mu)]), int(__LINE__, IK))

            sigmac = getMean(rand) - getOption(MU_DEF, mu)
            diff = abs(sigmac - getOption(SIGMA_DEF, sigma))
            reltol = getOption(SIGMA_DEF, sigma) * TOL
            assertion = assertion .and. diff < reltol
            call report(sigma, mu, sigmac)
            call test%assert(assertion, SK_"The condition `abs(sigmac - getOption(SIGMA_DEF, sigma)) < getOption(SIGMA_DEF, sigma) * TOL` must hold. Repeating the test may resolve the failure. present(sigma), present(mu) = "//getStr([present(sigma), present(mu)]), int(__LINE__, IK))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(sigma, mu, sigmac)
            integer :: i
            real(RKC), intent(in), optional :: sigma, mu, sigmac
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("present(mu)")
                call test%disp%show( present(mu) )
                if (present(mu)) then
                    call test%disp%show("mu")
                    call test%disp%show( mu )
                end if
                call test%disp%show("present(sigma)")
                call test%disp%show( present(sigma) )
                if (present(sigma)) then
                    call test%disp%show("sigma")
                    call test%disp%show( sigma )
                end if
                call test%disp%show("present(sigmac)")
                call test%disp%show( present(sigmac) )
                if (present(sigmac)) then
                    call test%disp%show("sigmac")
                    call test%disp%show( sigmac )
                end if
                call test%disp%show("reltol")
                call test%disp%show( reltol )
                call test%disp%show("diff")
                call test%disp%show( diff )
                do i = 1, NP
                    if (rand(i) < getOption(MU_DEF, mu)) then
                        call test%disp%show("[real(i, RKC), getOption(MU_DEF, mu), rand(i)]")
                        call test%disp%show( [real(i, RKC), getOption(MU_DEF, mu), rand(i)] )
                        exit
                    end if
                end do
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#else
            !%%%%%%%%%%%%%%%%%%%%%%%%
#error      "Unrecognized interface."
            !%%%%%%%%%%%%%%%%%%%%%%%%
#endif