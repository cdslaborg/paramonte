program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_parallelism, only: getImageID
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: getUnifRandState
    use pm_distUnif, only: rngf_type

    implicit none

    type(rngf_type) :: rng

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Set a random RNG seed.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("rng = rngf_type()")
                    rng = rngf_type()
    call disp%show("getUnifRandState()")
    call disp%show( getUnifRandState() )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Set the current RNG seed based on the input scalar `seed` and return a copy of seed vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("rng = rngf_type(seed = 12345_IK)")
                    rng = rngf_type(seed = 12345_IK)
    call disp%show("getUnifRandState()")
    call disp%show( getUnifRandState() )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Set the current RNG seed on each image distinctly based on the current date and time and and the image ID and return a copy of seed vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("rng = rngf_type(imageID = getImageID())")
                    rng = rngf_type(imageID = getImageID())
    call disp%show("getUnifRandState()")
    call disp%show( getUnifRandState() )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Set the current RNG seed on each image distinctly based on the input scalar `seed` and and the image ID and return a copy of seed vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("rng = rngf_type(seed = 12345_IK, imageID = getImageID())")
                    rng = rngf_type(seed = 12345_IK, imageID = getImageID())
    call disp%show("getUnifRandState()")
    call disp%show( getUnifRandState() )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Reset the RNG seed deterministically on each image distinctly based on the input scalar `seed` and and the image ID and return a copy of seed vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("rng = rngf_type(seed = 12345_IK, imageID = getImageID())")
                    rng = rngf_type(seed = 12345_IK, imageID = getImageID())
    call disp%show("getUnifRandState()")
    call disp%show( getUnifRandState() )
    call disp%show("getUnifRand(1, 100, 10)")
    call disp%show( getUnifRand(1, 100, 10) )
    call disp%show("rng = rngf_type(seed = 12345_IK, imageID = getImageID())")
                    rng = rngf_type(seed = 12345_IK, imageID = getImageID())
    call disp%show("getUnifRandState()")
    call disp%show( getUnifRandState() )
    call disp%show("getUnifRand(1, 100, 10)")
    call disp%show( getUnifRand(1, 100, 10) )
    call disp%show("rng = rngf_type() ! random RNG seed reset.")
                    rng = rngf_type() ! random RNG seed reset.
    call disp%show("getUnifRandState()")
    call disp%show( getUnifRandState() )
    call disp%show("getUnifRand(1, 100, 10)")
    call disp%show( getUnifRand(1, 100, 10) )
    call disp%skip()

end program example