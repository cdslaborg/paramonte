program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RK, RKS, RKD, RKH
    use pm_io, only: display_type
    use pm_distBern, only: getBernRand

    implicit none

    integer(IK) , parameter :: NP = 5_IK
    real(RKS)               :: rand_RKS(NP)
    real(RKD)               :: rand_RKD(NP)
    real(RKH)               :: rand_RKH(NP)
    type(display_type)      :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate `real`-valued random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar real random number.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("rand_RKS(1) = getBernRand(p = 0.5_RKS) ! 0.5 odds of success.")
                    rand_RKS(1) = getBernRand(p = 0.5_RKS)
    call disp%show("rand_RKS(1)")
    call disp%show( rand_RKS(1) )
    call disp%skip()

    call disp%show("rand_RKD(1) = getBernRand(p = 0.5_RKD)")
                    rand_RKD(1) = getBernRand(p = 0.5_RKD)
    call disp%show("rand_RKD(1)")
    call disp%show( rand_RKD(1) )
    call disp%skip()

    call disp%show("rand_RKH(1) = getBernRand(p = 0.5_RKH)")
                    rand_RKH(1) = getBernRand(p = 0.5_RKH)
    call disp%show("rand_RKH(1)")
    call disp%show( rand_RKH(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of real random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("rand_RKS(:) = getBernRand(p = 0.5_RKS, size = NP) ! 0.5 odds of success.")
                    rand_RKS(:) = getBernRand(p = 0.5_RKS, size = NP)
    call disp%show("rand_RKS(:)")
    call disp%show( rand_RKS(:) )
    call disp%skip()

    call disp%show("rand_RKD(:) = getBernRand(p = 0.5_RKD, size = NP)")
                    rand_RKD(:) = getBernRand(p = 0.5_RKD, size = NP)
    call disp%show("rand_RKD(:)")
    call disp%show( rand_RKD(:) )
    call disp%skip()

    call disp%show("rand_RKH(:) = getBernRand(p = 0.5_RKH, size = NP)")
                    rand_RKH(:) = getBernRand(p = 0.5_RKH, size = NP)
    call disp%show("rand_RKH(:)")
    call disp%show( rand_RKH(:) )
    call disp%skip()

end program example