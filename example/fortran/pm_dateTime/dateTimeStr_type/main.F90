program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_dateTime, only: dateTimeStr_type

    implicit none

    type(dateTimeStr_type) :: dts

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct a Gregorian calendar date and time.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("dts = dateTimeStr_type(year = '2020', zone = SK_'-5000', minute = '48')")
                    dts = dateTimeStr_type(year = '2020', zone = SK_'-5000', minute = '48')
    call disp%show("write(disp%unit, ""(*('''',g0,'''',:,', '))"") dts")
                    write(disp%unit, "(*('''',g0,'''',:,', '))") dts
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Construct the current date and time and moment in the Gregorian calendar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("dts = dateTimeStr_type()")
                    dts = dateTimeStr_type()
    call disp%show("write(disp%unit, ""(*('''',g0,'''',:,', '))"") dts")
                    write(disp%unit, "(*('''',g0,'''',:,', '))") dts
    call disp%skip()

end program example