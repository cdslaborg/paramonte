program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: getHour, getHour12

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("[getHour(), getHour12()]")
    call disp%show( [getHour(), getHour12()] )
    call disp%skip()

    call disp%skip()
    call disp%show("getHour12(0_IK)")
    call disp%show( getHour12(0_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHour12(1_IK)")
    call disp%show( getHour12(1_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHour12(11_IK)")
    call disp%show( getHour12(11_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHour12(12_IK)")
    call disp%show( getHour12(12_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHour12(13_IK)")
    call disp%show( getHour12(13_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHour12(23_IK)")
    call disp%show( getHour12(23_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("getHour12(24_IK)")
    call disp%show( getHour12(24_IK) )
    call disp%skip()

end program example