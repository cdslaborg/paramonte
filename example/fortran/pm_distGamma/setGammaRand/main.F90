program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKG => RKS ! all real kinds are supported.
    use pm_distGamma, only: setGammaRand
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    real(RKG), dimension(NP) :: Kappa, Sigma, rand

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLogSpace(Kappa, logx1 = log(0.1_RKG), logx2 = log(10._RKG))
    call setLogSpace(Sigma, logx1 = log(0.1_RKG), logx2 = log(10._RKG))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random numbers from the Gamma distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Gamma random value given integer shape and real inverse rate parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("Sigma(1)")
    call disp%show( Sigma(1) )
    call disp%show("call setGammaRand(rand(1:2), 1._RKG, sigma = Sigma(1))")
                    call setGammaRand(rand(1:2), 1._RKG, sigma = Sigma(1))
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Gamma random value given real shape and real inverse rate parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1)")
    call disp%show( Kappa(1) )
    call disp%show("Sigma(1)")
    call disp%show( Sigma(1) )
    call disp%show("call setGammaRand(rand(1:2), Kappa(1), sigma = Sigma(1))")
                    call setGammaRand(rand(1:2), Kappa(1), sigma = Sigma(1))
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Gamma random numbers with a fixed set of parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/3)")
    call disp%show( Kappa(1:NP:NP/3) )
    call disp%show("Sigma(1:NP:NP/3)")
    call disp%show( Sigma(1:NP:NP/3) )
    call disp%show("call setGammaRand(rand(1:NP:NP/3), Kappa(1), sigma = Sigma(1))")
                    call setGammaRand(rand(1:NP:NP/3), Kappa(1), sigma = Sigma(1))
    call disp%show("rand(1:NP:NP/3)")
    call disp%show( rand(1:NP:NP/3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Gamma random numbers for a range of parameters.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Kappa(1:NP:NP/3)")
    call disp%show( Kappa(1:NP:NP/3) )
    call disp%show("Sigma(1:NP:NP/3)")
    call disp%show( Sigma(1:NP:NP/3) )
    call disp%show("call setGammaRand(rand(1:NP:NP/3), Kappa(1:NP:NP/3), sigma = Sigma(1:NP:NP/3))")
                    call setGammaRand(rand(1:NP:NP/3), Kappa(1:NP:NP/3), sigma = Sigma(1:NP:NP/3))
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
        do itry = 1, 30
            call disp%skip()
            call disp%show("kappa = getExpRand(1._RKG); sigma = getExpRand(1._RKG)")
                            kappa = getExpRand(1._RKG); sigma = getExpRand(1._RKG)
            call disp%show("[kappa, sigma]")
            call disp%show( [kappa, sigma] )
            call disp%show("call setGammaRand(rand, kappa, sigma)")
                            call setGammaRand(rand, kappa, sigma)
            call disp%show("mean = kappa * sigma")
                            mean = kappa * sigma
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
        real(RKG):: rand(5000, 3)
        call setGammaRand(rand(:, 1), +0.8_RKG, sigma = 2._RKG)
        call setGammaRand(rand(:, 2), +1.0_RKG, sigma = 2._RKG)
        call setGammaRand(rand(:, 3), +5.0_RKG, sigma = 2._RKG)
        if (0 /= getErrTableWrite(SK_"setGammaRand.RK.txt", rand)) error stop "Table writing failed."
    end block

end program example