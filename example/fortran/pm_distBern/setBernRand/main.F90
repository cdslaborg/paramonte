program example

    use pm_kind, only: SK, IK, LK ! all intrinsic types and kinds are supported.
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_distBern, only: setBernRand

    implicit none

    integer(IK) , parameter :: NP = 5_IK
    logical                 :: Rand_LK(NP)
    integer                 :: Rand_IK(NP)
    real                    :: Rand_RK(NP)
    type(display_type)      :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate `integer`-valued random state.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar integer random value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setBernRand(Rand_IK(1), getUnifRand(0., 1.), p = .5) ! 0.5 odds of success.")
                    call setBernRand(Rand_IK(1), getUnifRand(0., 1.), p = .5)
    call disp%show("Rand_IK(1)")
    call disp%show( Rand_IK(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector integer random value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setBernRand(Rand_IK(:), getUnifRand(0., 1., NP), p = 0.75) ! 0.75 odds of success.")
                    call setBernRand(Rand_IK(:), getUnifRand(0., 1., NP), p = 0.75)
    call disp%show("Rand_IK(:)")
    call disp%show( Rand_IK(:) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate `logical`-valued random state.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar logical random value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setBernRand(Rand_LK(1), getUnifRand(0., 1.), p = .5) ! 0.5 odds of success.")
                    call setBernRand(Rand_LK(1), getUnifRand(0., 1.), p = .5)
    call disp%show("Rand_LK(1)")
    call disp%show( Rand_LK(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector logical random value.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setBernRand(Rand_LK(:), getUnifRand(0., 1., NP), p = 0.75) ! 0.75 odds of success.")
                    call setBernRand(Rand_LK(:), getUnifRand(0., 1., NP), p = 0.75)
    call disp%show("Rand_LK(:)")
    call disp%show( Rand_LK(:) )
    call disp%skip()

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

    call disp%show("call setBernRand(Rand_RK(1), getUnifRand(0., 1.), p = .5) ! 0.5 odds of success.")
                    call setBernRand(Rand_RK(1), getUnifRand(0., 1.), p = .5)
    call disp%show("Rand_RK(1)")
    call disp%show( Rand_RK(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of real random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setBernRand(Rand_RK(:), getUnifRand(0., 1., NP), p = 0.75) ! 0.75 odds of success.")
                    call setBernRand(Rand_RK(:), getUnifRand(0., 1., NP), p = 0.75)
    call disp%show("Rand_RK(:)")
    call disp%show( Rand_RK(:) )
    call disp%skip()

end program example