program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRandStateSize

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Get the RNG seed size.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getUnifRandStateSize()")
    call disp%show( getUnifRandStateSize() )
    call disp%skip()

end program example