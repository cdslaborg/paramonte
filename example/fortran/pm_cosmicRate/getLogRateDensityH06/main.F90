program example

    use pm_kind, only: SK, IK, LK, RKG => RKD
    use pm_cosmicRate, only: getLogRateDensityH06
    use pm_arraySpace, only: getLinSpace
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    real(RKG) :: redshift = 5.5_RKG
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the log(RateDensity) according to the H06 LGRB rate density parameters for the Hopkins and Beacom (2006) SFR model.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("redshift")
    call disp%show( redshift )
    call disp%show("getLogRateDensityH06(log(redshift + 1))")
    call disp%show( getLogRateDensityH06(log(redshift + 1)) )
    call disp%skip()

    ! Generate both the cosmic rate and the rate density.

    block
        use pm_cosmology, only: getVolComDiffNormed
        integer(IK) :: fileUnit, i
        real(RKG) :: maxRateFormLGRB, maxRateDensityFormLGRB
        real(RKG), allocatable :: zplus1(:), logzplus1(:), rateFormLGRB(:), rateDensityFormLGRB(:)
        logzplus1 = getLinSpace(0.01_RKG, log(13._RKG), 500_IK)
        zplus1 = exp(logzplus1)
        rateDensityFormLGRB = exp(getLogRateDensityH06(logzplus1))
        rateFormLGRB = rateDensityFormLGRB * getVolComDiffNormed(zplus1, reltol = sqrt(epsilon(0._RKG)))
        maxRateDensityFormLGRB = maxval(rateDensityFormLGRB)
        maxRateFormLGRB = maxval(rateFormLGRB)
        open(newunit = fileUnit, file = "getLogRateDensityH06.csv")
        write(fileUnit, "(*(g0,:,','))") "redshift, rateFormLGRB, rateDensityFormLGRB"
        do i = 1, size(zplus1)
            write(fileUnit, "(*(g0,:,','))") zplus1(i) - 1, rateFormLGRB(i) / maxRateFormLGRB, rateDensityFormLGRB(i) / maxRateDensityFormLGRB
        end do
        close(fileUnit)
    end block

    ! Set up a Markov Chain Monte Carlo sampler to generate random sample from the target rate density model.

    block
        use pm_err, only: err_type
        use pm_sampling, only: paradram_type, getErrSampling
        type(paradram_type) :: sampler
        type(err_type) :: err
        sampler%outputFileName = "./zdistH06"
        sampler%outputStatus = "retry"
        sampler%domainAxisName = ["redshift"]
        sampler%domainCubeLimitLower = [0._RKG]
        sampler%outputSampleSize = 2500
        sampler%outputChainSize = 5000
        sampler%proposalStart = [3]
        err = getErrSampling(sampler, getLogFunc, 1_IK)
        if (err%occurred) error stop err%msg
    end block

contains

    recursive function getLogFunc(redshift) result(logRateDensity)
        use pm_cosmology, only: getVolComDiffNormed
        real(RKG), intent(in), contiguous :: redshift(:)
        real(RKG) :: logRateDensity
        real(RKG) :: zplus1
        zplus1 = redshift(1) + 1
        logRateDensity = getLogRateDensityH06(log(zplus1)) + log(getVolComDiffNormed(zplus1, reltol = sqrt(epsilon(0._RKG))))
    end function

end program example