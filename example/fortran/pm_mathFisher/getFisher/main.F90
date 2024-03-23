program example

    use pm_kind, only: SK, IK, LK, RKC => RKD
    use pm_arraySpace, only: getLinSpace
    use pm_mathFisher, only: getFisher
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    real(RKC) :: redshift = 5.5_RKC
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Fisher z transformation of a random bounded variable.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKC => RKS
        real(RKC), allocatable :: rand(:)
        call disp%skip()
        call disp%show("rand = getUnifRand(0._RKC, 1._RKC, 5_IK)")
                        rand = getUnifRand(0._RKC, 1._RKC, 5_IK)
        call disp%show("rand")
        call disp%show( rand )
        call disp%show("getFisher(rand)")
        call disp%show( getFisher(rand) )
        call disp%show("getFisher(rand, lb = -1._RKC, ub = 1._RKC)")
        call disp%show( getFisher(rand, lb = -1._RKC, ub = 1._RKC) )
        call disp%show("getFisher(rand, lb = 0._RKC, ub = 1._RKC)")
        call disp%show( getFisher(rand, lb = 0._RKC, ub = 1._RKC) )
        call disp%skip()
    end block

    ! Generate both the cosmic rate and the rate density.

    block
        use pm_io, only: getErrTableWrite
        real(RKC), allocatable :: rho(:)
        rho = getLinSpace(-.99, +.99, 500_IK)
        if (0 /= getErrTableWrite("getFisher.csv", reshape([rho, getFisher(rho)], [size(rho), 2]), header = "Correlation Coefficient,Fisher Transformation")) error stop "Table writing failed."
    end block

end program example