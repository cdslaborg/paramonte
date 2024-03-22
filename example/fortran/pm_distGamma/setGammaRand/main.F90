program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RK => RKS ! all real kinds are supported.
    use pm_distGamma, only: setGammaRand
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    real(RK), dimension(NP) :: Kappa, Sigma, rand

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLogSpace(Kappa, logx1 = log(0.1_RK), logx2 = log(10._RK))
    call setLogSpace(Sigma, logx1 = log(0.1_RK), logx2 = log(10._RK))

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
    call disp%show("call setGammaRand(rand(1:2), 1._RK, sigma = Sigma(1))")
                    call setGammaRand(rand(1:2), 1._RK, sigma = Sigma(1))
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

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        integer(IK), parameter :: NP = 5000_IK
        real(RK), dimension(NP) :: Rand1, Rand2, Rand3
        call setGammaRand(Rand1, +0.8_RK, sigma = 2._RK)
        call setGammaRand(Rand2, +1.0_RK, sigma = 2._RK)
        call setGammaRand(Rand3, +5.0_RK, sigma = 2._RK)
        open(newunit = fileUnit, file = "setGammaRand.RK.txt")
        write(fileUnit,"(3(g0,:,' '))") ( Rand1(i) &
                                        , Rand2(i) &
                                        , Rand3(i) &
                                        , i = 1,NP &
                                        )
        close(fileUnit)
    end block

end program example