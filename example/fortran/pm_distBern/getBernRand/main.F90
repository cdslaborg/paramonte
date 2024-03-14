program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RK, RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_distBern, only: getBernRand

    implicit none

    integer(IK) , parameter :: NP = 5_IK
    real(RK32 )             :: Rand_RK32 (NP)
    real(RK64 )             :: Rand_RK64 (NP)
    real(RK128)             :: Rand_RK128(NP)
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

    call disp%show("rand_RK32(1) = getBernRand(p = 0.5_RK32) ! 0.5 odds of success.")
                    rand_RK32(1) = getBernRand(p = 0.5_RK32)
    call disp%show("rand_RK32(1)")
    call disp%show( rand_RK32(1) )
    call disp%skip()

    call disp%show("rand_RK64(1) = getBernRand(p = 0.5_RK64)")
                    rand_RK64(1) = getBernRand(p = 0.5_RK64)
    call disp%show("rand_RK64(1)")
    call disp%show( rand_RK64(1) )
    call disp%skip()

    call disp%show("rand_RK128(1) = getBernRand(p = 0.5_RK128)")
                    rand_RK128(1) = getBernRand(p = 0.5_RK128)
    call disp%show("rand_RK128(1)")
    call disp%show( rand_RK128(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of real random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("rand_RK32(:) = getBernRand(p = 0.5_RK32, size = NP) ! 0.5 odds of success.")
                    rand_RK32(:) = getBernRand(p = 0.5_RK32, size = NP)
    call disp%show("rand_RK32(:)")
    call disp%show( rand_RK32(:) )
    call disp%skip()

    call disp%show("rand_RK64(:) = getBernRand(p = 0.5_RK64, size = NP)")
                    rand_RK64(:) = getBernRand(p = 0.5_RK64, size = NP)
    call disp%show("rand_RK64(:)")
    call disp%show( rand_RK64(:) )
    call disp%skip()

    call disp%show("rand_RK128(:) = getBernRand(p = 0.5_RK128, size = NP)")
                    rand_RK128(:) = getBernRand(p = 0.5_RK128, size = NP)
    call disp%show("rand_RK128(:)")
    call disp%show( rand_RK128(:) )
    call disp%skip()

end program example