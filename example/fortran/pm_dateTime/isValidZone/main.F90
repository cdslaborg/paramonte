program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: isValidZone
    use pm_val2str, only: getStr

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isValidZone(0_IK)")
    call disp%show( isValidZone(0_IK) )
    call disp%skip()

    call disp%skip()
    call disp%show("isValidZone(60_IK * [integer(IK) :: -13, -12, 12, 14, 16])")
    call disp%show( isValidZone(60_IK * [integer(IK) :: -13, -12, 12, 14, 16]) )
    call disp%skip()

end program example