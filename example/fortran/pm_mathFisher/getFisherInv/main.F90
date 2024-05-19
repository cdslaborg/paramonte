program example

    use pm_kind, only: SK, IK, LK, RKG => RKD
    use pm_arraySpace, only: getLinSpace
    use pm_mathFisher, only: getFisherInv
    use pm_mathFisher, only: getFisher
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    real(RKG) :: redshift = 5.5_RKG
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Inverse Fisher z transformation of a random bounded variable.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKS
        real(RKG), allocatable :: rand(:)
        call disp%skip()
        call disp%show("rand = getUnifRand(-9, 9, 5_IK)")
                        rand = getUnifRand(-9, 9, 5_IK)
        call disp%show("rand")
        call disp%show( rand )
        call disp%show("getFisherInv(rand)")
        call disp%show( getFisherInv(rand) )
        call disp%show("getFisherInv(rand, lb = -1._RKG, ub = 1._RKG)")
        call disp%show( getFisherInv(rand, lb = -1._RKG, ub = 1._RKG) )
        call disp%show("getFisherInv(rand, lb = 0._RKG, ub = 1._RKG)")
        call disp%show( getFisherInv(rand, lb = 0._RKG, ub = 1._RKG) )
        call disp%skip()
    end block

    ! Generate both the cosmic rate and the rate density.

    block
        use pm_io, only: getErrTableWrite
        real(RKG), allocatable :: fish(:)
        fish = getFisher(getLinSpace(-.99, +.99, 500_IK))
        if (0 /= getErrTableWrite("getFisherInv.csv", reshape([fish, getFisherInv(fish)], [size(fish), 2]), header = "Fisher Transformation,Correlation Coefficient")) error stop "Table writing failed."
    end block

end program example