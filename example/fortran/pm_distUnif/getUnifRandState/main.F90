program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_parallelism, only: getImageID
    use pm_distUnif, only: getUnifRandState

    implicit none

    integer(IK), allocatable :: unifRandState(:)

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the current RNG seed.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("unifRandState = getUnifRandState()")
                    unifRandState = getUnifRandState()
    call disp%show("unifRandState")
    call disp%show( unifRandState )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Set the current RNG seed based on the input scalar `seed` and return a copy of seed vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("unifRandState = getUnifRandState(seed = 12345_IK)")
                    unifRandState = getUnifRandState(seed = 12345_IK)
    call disp%show("unifRandState")
    call disp%show( unifRandState )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Set the current RNG seed on each image distinctly based on the current date and time and and the image ID and return a copy of seed vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("unifRandState = getUnifRandState(imageID = getImageID())")
                    unifRandState = getUnifRandState(imageID = getImageID())
    call disp%show("unifRandState")
    call disp%show( unifRandState )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Set the current RNG seed on each image distinctly based on the input scalar `seed` and and the image ID and return a copy of seed vector.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("unifRandState = getUnifRandState(seed = 12345_IK, imageID = getImageID())")
                    unifRandState = getUnifRandState(seed = 12345_IK, imageID = getImageID())
    call disp%show("unifRandState")
    call disp%show( unifRandState )
    call disp%skip()

end program example