program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getSecond

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getSecond()")
    call disp%show( getSecond() )
    call disp%skip()

end program example