program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RKS ! all real kinds are supported.
    use pm_distGenGamma, only: setGenGammaRand
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 10000_IK
    real(RKG), dimension(NP) :: kappa, omega, sigma, rand, mean

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLogSpace(kappa, logx1 = log(0.1_RKG), logx2 = log(10._RKG))
    call setLogSpace(sigma, logx1 = log(0.1_RKG), logx2 = log(10._RKG))
    call setLogSpace(omega, logx1 = log(0.1_RKG), logx2 = log(10._RKG))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random numbers from the GenGamma distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("sigma(1:NP:NP/3)")
    call disp%show( sigma(1:NP:NP/3) )
    call disp%show("call setGenGammaRand(rand(1:NP:NP/3), kappa(1), omega(1), sigma(1:NP:NP/3))")
                    call setGenGammaRand(rand(1:NP:NP/3), kappa(1), omega(1), sigma(1:NP:NP/3))
    call disp%show("rand(1:NP:NP/3)")
    call disp%show( rand(1:NP:NP/3) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1:NP:NP/3)")
    call disp%show( kappa(1:NP:NP/3) )
    call disp%show("call setGenGammaRand(rand(1:NP:NP/3), kappa(1:NP:NP/3), omega(1), sigma(1))")
                    call setGenGammaRand(rand(1:NP:NP/3), kappa(1:NP:NP/3), omega(1), sigma(1))
    call disp%show("rand(1:NP:NP/3)")
    call disp%show( rand(1:NP:NP/3) )
    call disp%skip()

    call disp%skip()
    call disp%show("kappa(1:NP:NP/3)")
    call disp%show( kappa(1:NP:NP/3) )
    call disp%show("sigma(1:NP:NP/3)")
    call disp%show( sigma(1:NP:NP/3) )
    call disp%show("call setGenGammaRand(rand(1:NP:NP/3), kappa(1:NP:NP/3), omega(1), sigma(1:NP:NP/3))")
                    call setGenGammaRand(rand(1:NP:NP/3), kappa(1:NP:NP/3), omega(1), sigma(1:NP:NP/3))
    call disp%show("rand(1:NP:NP/3)")
    call disp%show( rand(1:NP:NP/3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Test the mean of a random sample against the analytic answer.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_sampleMean, only: getMean
        use pm_distExp, only: getExpRand
        real(RKG) :: kappa, omega, sigma, mean
        integer(IK) :: itry
        do itry = 1, 10
            call disp%skip()
            call disp%show("kappa = getExpRand(1._RKG); omega = getExpRand(1._RKG); sigma = getExpRand(1._RKG)")
                            kappa = getExpRand(1._RKG); omega = getExpRand(1._RKG); sigma = getExpRand(1._RKG)
            call disp%show("[kappa, omega, sigma]")
            call disp%show( [kappa, omega, sigma] )
            call disp%show("call setGenGammaRand(rand, kappa, omega, sigma)")
                            call setGenGammaRand(rand, kappa, omega, sigma)
            call disp%show("mean = exp(log_gamma(kappa + omega) - log_gamma(kappa)) * sigma")
                            mean = exp(log_gamma(kappa + omega) - log_gamma(kappa)) * sigma
            call disp%show("[getMean(rand), mean]")
            call disp%show( [getMean(rand), mean] )
            call disp%skip()
        end do
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_io, only: getErrTableWrite
        integer(IK) :: fileUnit, i
        integer(IK), parameter :: NP = 5000
        real(RKG), dimension(NP, 5) :: rand
        call setGenGammaRand(rand(:,1), 2.0_RKG, omega = 1 / 5.0_RKG, sigma = 1 / 0.3_RKG)
        call setGenGammaRand(rand(:,2), .14_RKG, omega = 1 / 7.0_RKG, sigma = 1 / .14_RKG)
        call setGenGammaRand(rand(:,3), 0.2_RKG, omega = 1 / 5.0_RKG, sigma = 1 / 0.2_RKG)
        call setGenGammaRand(rand(:,4), 0.5_RKG, omega = 1 / 2.0_RKG, sigma = 1 / 0.5_RKG)
        call setGenGammaRand(rand(:,5), 2.0_RKG, omega = 1 / 0.5_RKG, sigma = 1 / 1.0_RKG)
        !call setGenGammaRand(rand(:,6), 1.0_RKG, omega = 1 / 0.5_RKG, sigma = 1 / 0.5_RKG)
        if (0 /= getErrTableWrite("setGenGammaRand.RKG.txt", rand)) error stop 'table write failed.'
    end block

end program example